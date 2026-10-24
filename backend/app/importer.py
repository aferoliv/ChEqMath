"""Importa componentes_e_especies_ChEq_Math.xlsx para SQLite (executar uma vez por versão da planilha).

Regras, de acordo com Teoria.docx:
- aba 'components': nome (A) e código (C). 0001 = água (solvente).
- aba 'Species': espécies explícitas. Colunas: nome (A), tipo (B), componentes (C:E),
  coeficientes ν (G:I), log K original (K), log K (L). Usa-se a coluna 'log K' (L).
- aba 'pKa': para cada ligante, geram-se as espécies protonadas H_k·X = X + k·H+,
  com log β_k = soma dos k maiores pKa (valores em ordem crescente). Sistemas
  metálicos (nome com carga, ex.: 'Cd(II)') não são gerados: a hidrólise vem da aba Species.
- Espécies com mesmos componentes e coeficientes são deduplicadas (prevalece a explícita).
- Consistência de sinal: se a espécie contém água e H+ com ν < 0 (hidrólise), log K deve ser
  negativo. O hidróxido (água +1, H+ −1) com +15,75 é corrigido para −pKw (25 °C, I = 0).
  Outros casos inconsistentes são apenas registrados como alerta.

Uso:
    python -m app.importer [--xlsx CAMINHO] [--db CAMINHO]
"""
from __future__ import annotations

import argparse
import os
import re
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path

import openpyxl

from .engine import PROTON, SOLVENT

# Regra: toda constante de espécie que contém água (ν_água ≠ 0) está referida à água em 55,6 mol L⁻¹
# (ver engine.effective_log_beta). Constantes sem água não são alteradas.

METAL_RE = re.compile(r"\((I{1,3}|\d\+)\)")
SKIP_SYSTEMS = {"acid - strong", "base - strong", "test"}
ALIASES = {"acetic acid": "acetate"}  # nome na aba pKa -> nome do componente
SCHEMA = Path(__file__).with_name("schema.sql")
REPO = Path(__file__).resolve().parents[2]
DEFAULT_XLSX = REPO / "data" / "componentes_e_especies_ChEq_Math.xlsx"
DEFAULT_DB = REPO / "data" / "speciation.sqlite"


def _norm(text) -> str:
    return re.sub(r"\s+", " ", str(text or "")).strip().lower()


def _num(value):
    """Converte número da planilha; aceita vírgula decimal ('-6,7')."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip().replace(",", ".")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _sheet(wb, prefix: str):
    for ws in wb.worksheets:
        if ws.title.strip().lower().startswith(prefix):
            return ws
    raise KeyError(f"aba '{prefix}' não encontrada na planilha")


def _sign_ok(nu: dict, log_beta: float) -> bool:
    h = nu.get(PROTON, 0)
    if h == 0 or nu.get(SOLVENT, 0) == 0:
        return True
    return (log_beta < 0) == (h < 0)


@dataclass
class Built:
    components: list = field(default_factory=list)  # dicts: code, name, smiles, is_solvent
    species: list = field(default_factory=list)     # dicts: code, name, kind, log_beta, log_beta_raw, nu, cap_of, source, note
    pka: list = field(default_factory=list)         # dicts: comp_code, idx, value, carga_max
    log: list = field(default_factory=list)         # (nível, mensagem)


def _read_components(ws) -> list:
    out = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        name, smiles, code = row[0], row[1], row[2]
        if not name or not code:
            continue
        code = str(code).strip()
        out.append({
            "code": code,
            "name": str(name).strip(),
            "smiles": str(smiles).strip() if smiles else None,
            "is_solvent": int(code == SOLVENT),
        })
    return out


def _read_explicit_species(ws, b: Built) -> list:
    recs = []
    for r_i, row in enumerate(ws.iter_rows(min_row=2, values_only=True), start=2):
        name, kind = row[0], row[1]
        if not name or not kind:
            continue  # cabeçalho ou linha de rótulo 'species'
        nu, ok = [], True
        for comp, coef in zip(row[2:5], row[6:9]):
            if comp is None or str(comp).strip() == "":
                continue
            value = _num(coef)
            if value is None:
                b.log.append(("erro", f"Species linha {r_i}: coeficiente ausente para {comp}"))
                ok = False
                break
            nu.append((str(comp).strip(), int(value)))
        raw = _num(row[11])
        if not ok:
            continue
        if raw is None:
            b.log.append(("erro", f"Species linha {r_i}: log K inválido ({row[11]!r})"))
            continue
        str_kind = str(kind).strip()
        name = str(name).strip()
        nu_t = tuple(nu)
        # Constantes são usadas como estão na coluna log K. Erros são corrigidos na planilha.
        if not _sign_ok(dict(nu_t), raw):
            b.log.append(("alerta", f"{name} (linha {r_i}): log K = {raw:+.2f} com sinal inconsistente "
                                    f"com a estequiometria; verificar na planilha"))
        recs.append({"name": name, "kind": str_kind, "nu": nu_t, "log_beta": raw,
                     "log_beta_raw": raw, "source": "planilha:Species", "note": None, "row": r_i})
    return recs


def _hydrolysis_series_alert(recs: list, b: Built) -> None:
    """Alerta (sem alterar valores) para séries M(OH)n com n ≥ 2.

    Os valores são usados como estão. Se forem pK sucessivos (stepwise) e não constantes
    cumulativas, a correção deve ser feita na planilha.
    """
    groups: dict = {}
    for r in recs:
        if r["kind"] != "ABB":
            continue
        d = dict(r["nu"])
        metals = [c for c, v in d.items() if c not in (SOLVENT, PROTON) and v > 0]
        n = -d.get(PROTON, 0)
        if len(metals) == 1 and n >= 1 and d.get(SOLVENT, 0) == n and d.get(metals[0]) == 1:
            groups.setdefault(metals[0], []).append((n, r["name"]))
    for metal, items in groups.items():
        if len(items) >= 2:
            names = ", ".join(f"{nm} (n={n})" for n, nm in sorted(items))
            b.log.append(("alerta", f"série de hidrólise {names}: confirmar se log K são cumulativos "
                                    f"(log β_n) ou stepwise (pK sucessivos)"))


def _generate_from_pka(ws, comp_by_name: dict, existing_nu: set, b: Built) -> list:
    recs = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        name = row[0]
        if not name:
            continue
        key = _norm(name)
        carga = _num(row[1])
        vals = [v for v in (_num(x) for x in row[2:]) if v is not None]
        if key in SKIP_SYSTEMS:
            b.log.append(("info", f"pKa '{str(name).strip()}': não é componente (ignorado)"))
            continue
        if METAL_RE.search(str(name)):
            b.log.append(("info", f"pKa '{str(name).strip()}': hidrólise de metal, use a aba Species"))
            continue
        code = comp_by_name.get(ALIASES.get(key, key))
        if code is None:
            b.log.append(("alerta", f"pKa '{str(name).strip()}': sem componente correspondente"))
            continue
        if not vals:
            continue
        vals = sorted(vals)
        n = len(vals)
        for idx, v in enumerate(vals, start=1):
            b.pka.append({"comp_code": code, "idx": idx, "value": v,
                          "carga_max": int(carga) if carga is not None else None})
        for k in range(1, n + 1):
            nu = ((code, 1), (PROTON, k))
            key_nu = frozenset(dict(nu).items())
            if key_nu in existing_nu:
                continue
            existing_nu.add(key_nu)
            log_beta = float(sum(vals[n - k:]))
            recs.append({"name": f"{str(name).strip()} H{k}", "kind": "ABB", "nu": nu,
                         "log_beta": log_beta, "log_beta_raw": log_beta,
                         "source": "gerada:pKa", "note": f"H{k}·X de pKa {', '.join(f'{v:.3f}' for v in vals)}",
                         "row": None})
    return recs


def build(xlsx_path: Path = DEFAULT_XLSX) -> Built:
    b = Built()
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    b.components = _read_components(_sheet(wb, "components"))
    comp_by_name: dict = {}
    for c in b.components:
        comp_by_name.setdefault(_norm(c["name"]), c["code"])

    explicit = _read_explicit_species(_sheet(wb, "species"), b)
    _hydrolysis_series_alert(explicit, b)
    aq_recs = [r for r in explicit if r["kind"] != "S0"]
    cap_recs = [r for r in explicit if r["kind"] == "S0"]
    existing_nu = {frozenset(dict(r["nu"]).items()) for r in aq_recs}
    aq_recs += _generate_from_pka(_sheet(wb, "pka"), comp_by_name, existing_nu, b)

    for i, r in enumerate(aq_recs, start=1):
        r["code"] = f"SP{i:04d}"
        r["cap_of"] = None
        r["water_ref"] = 1 if SOLVENT in dict(r["nu"]) else 0
        b.species.append(r)

    next_i = len(aq_recs) + 1
    for r in cap_recs:
        target = next(
            (a for a in aq_recs
             if _norm(a["name"]) == _norm(r["name"]) and dict(a["nu"]) == dict(r["nu"])),
            None,
        )
        if target is None:
            b.log.append(("erro", f"S0 '{r['name']}' sem espécie aquosa de mesmos componentes: ignorado"))
            continue
        r["code"] = f"SP{next_i:04d}"
        next_i += 1
        r["cap_of"] = target["code"]
        r["log_beta"] = r["log_beta"]  # log So
        r["water_ref"] = 0
        b.species.append(r)

    return b


def write_db(b: Built, db_path: Path) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    try:
        conn.executescript(SCHEMA.read_text(encoding="utf-8"))
        conn.executemany(
            "INSERT INTO components(code, name, smiles, is_solvent) VALUES (?, ?, ?, ?)",
            [(c["code"], c["name"], c["smiles"], c["is_solvent"]) for c in b.components],
        )
        conn.executemany(
            "INSERT INTO species(code, name, kind, log_beta, log_beta_raw, cap_of, water_ref, source, note) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [(s["code"], s["name"], s["kind"], s["log_beta"], s["log_beta_raw"],
              s["cap_of"], s.get("water_ref", 0), s["source"], s["note"]) for s in b.species],
        )
        conn.executemany(
            "INSERT INTO stoich(species_code, comp_code, nu) VALUES (?, ?, ?)",
            [(s["code"], j, nu) for s in b.species for j, nu in s["nu"]],
        )
        conn.executemany(
            "INSERT INTO pka(comp_code, idx, value, carga_max) VALUES (?, ?, ?, ?)",
            [(p["comp_code"], p["idx"], p["value"], p["carga_max"]) for p in b.pka],
        )
        conn.executemany("INSERT INTO import_log(level, message) VALUES (?, ?)", b.log)
        conn.commit()
    finally:
        conn.close()


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Importa a planilha ChEq para SQLite")
    p.add_argument("--xlsx", default=os.environ.get("XLSX_PATH", str(DEFAULT_XLSX)))
    p.add_argument("--db", default=os.environ.get("DB_PATH", str(DEFAULT_DB)))
    args = p.parse_args(argv)
    b = build(Path(args.xlsx))
    write_db(b, Path(args.db))
    n_alert = sum(1 for lvl, _ in b.log if lvl in ("alerta", "erro"))
    print(f"componentes: {len(b.components)} | espécies: {len(b.species)} | pKa: {len(b.pka)} | alertas/erros: {n_alert}")
    for level, msg in b.log:
        print(f"[{level}] {msg}")


if __name__ == "__main__":
    main()
