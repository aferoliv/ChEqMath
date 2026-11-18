"""Importa componentes_e_especies_ChEq_Math.xlsx para SQLite (executar uma vez por versão da planilha).

Regras, de acordo com Teoria.docx:
- aba 'components': nome (A) e código (C). 0001 = água (solvente).
- aba 'Species': espécies explícitas. Colunas: nome (A), tipo (B), componentes (C:E),
  coeficientes ν (G:I), log K original (K), log K (L). Usa-se a coluna 'log K' (L).
- aba 'pKa': NÃO é usada. Serve apenas de apoio ao usuário; todas as espécies vêm da aba Species.
- Espécies com mesmos componentes e coeficientes são deduplicadas (prevalece a explícita).
- Consistência de sinal: casos inconsistentes são apenas registrados como alerta.
- Carga dos componentes (z_j, para Debye–Hückel). A coluna de carga da aba components
  (cabeçalho 'Maximum Electrical Charge', ou contendo 'charge'/'carga') traz a
  CARGA MÁXIMA DO SISTEMA associado ao componente, isto é, a carga da espécie totalmente
  protonada (Cd(II): '2+'; glicina: '+', de H2Gly+; acetato: '0', de HAc). Aceita '2+',
  '+', '-', '3-', '+2', '-1', '0' ou número. Como o componente é a forma totalmente
  desprotonada:
      z_j = q_max,j − n_H,j
  com n_H,j = número de prótons ácidos do sistema = maior valor entre (a) o nº de pKa da aba
  pKa (usada só para isso) e (b) o maior k das espécies de protonação pura HₖX
  (ν = {X: 1, H+: k}) da aba Species. Metais: n_H = 0 (hidrólise envolve água), z = q_max.
  Água: 0; H+: +1.
  As colunas da aba components são localizadas pelo cabeçalho ('code', 'category …').
  Célula de carga vazia: metal pelo nome ('(II)', '(2+)') com alerta; sem isso, carga desconhecida
  (NULL) e alerta — o cálculo com Debye–Hückel recusa esse componente.

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
ROMAN = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}
CHARGE_ROMAN_RE = re.compile(r"\(\s*(VI|IV|V|I{1,3})\s*\)")
CHARGE_PLUS_RE = re.compile(r"\(\s*(\d)\s*\+\s*\)")
CHARGE_TXT_RE = re.compile(r"^([+-]?)(\d*)([+-]?)$")
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
    components: list = field(default_factory=list)  # dicts: code, name, is_solvent, charge_max, charge, categories
    species: list = field(default_factory=list)     # dicts: code, name, kind, log_beta, log_beta_raw, nu, cap_of, source, note
    log: list = field(default_factory=list)         # (nível, mensagem)


def _parse_charge(value):
    """'2+' → 2, '+' → 1, '-' → −1, '3-' → −3, '+2' → 2, '-1' → −1, '0' → 0, 2.0 → 2.

    Vazio → None. Texto inválido → ValueError.
    """
    if value is None:
        return None
    if isinstance(value, (int, float)):
        if float(value) != int(value):
            raise ValueError(f"carga não inteira: {value!r}")
        return int(value)
    text = str(value).strip().replace("−", "-").replace(" ", "")
    if not text:
        return None
    m = CHARGE_TXT_RE.match(text)
    if not m or (m.group(1) and m.group(3)):
        raise ValueError(f"carga inválida: {value!r}")
    pre, digits, post = m.groups()
    sign = pre or post
    if not digits:
        if not sign:
            raise ValueError(f"carga inválida: {value!r}")
        return 1 if sign == "+" else -1
    n = int(digits)
    return -n if sign == "-" else n


def _component_columns(header: list) -> dict:
    """Localiza as colunas da aba components pelo cabeçalho (linha 1).

    nome: coluna A; carga: cabeçalho com 'charge' ou 'carga'; código: 'code'/'código';
    categorias: cabeçalhos que começam com 'categor'. Sem cabeçalho reconhecível, usa o
    leiaute antigo (A nome, B carga, C código, D:F categorias).
    """
    h = [_norm(x) for x in header]
    def find(pred):
        return next((k for k, x in enumerate(h) if x and pred(x)), None)
    code = find(lambda x: x in ("code", "código", "codigo"))
    charge = find(lambda x: "charge" in x or "carga" in x)
    cats = [k for k, x in enumerate(h) if x.startswith("categor")]
    if code is None:
        return {"name": 0, "charge": 1, "code": 2, "cats": [3, 4, 5]}
    return {"name": 0, "charge": charge, "code": code, "cats": cats}


def _read_components(ws, b: "Built") -> list:
    """Aba components: nome, carga máxima do sistema, código e categorias (por cabeçalho)."""
    header = next(ws.iter_rows(min_row=1, max_row=1, values_only=True))
    col = _component_columns(list(header))
    if col["charge"] is None:
        b.log.append(("alerta", "aba components sem coluna de carga máxima ('Maximum Electrical "
                                "Charge'/'carga'): cargas só por reserva"))
    get = lambda row, k: row[k] if k is not None and k < len(row) else None
    out = []
    for r_i, row in enumerate(ws.iter_rows(min_row=2, values_only=True), start=2):
        name, qcell, code = get(row, col["name"]), get(row, col["charge"]), get(row, col["code"])
        if not name or not code:
            continue
        code = str(code).strip()
        try:
            qmax = _parse_charge(qcell)
        except ValueError as exc:
            b.log.append(("erro", f"components linha {r_i} ({str(name).strip()}): {exc}; ignorada"))
            qmax = None
        cats = []
        for k in col["cats"]:
            cell = get(row, k)
            if cell is None or str(cell).strip() == "":
                continue
            cat = str(cell).strip().lower()
            if cat not in cats:
                cats.append(cat)
        out.append({
            "code": code,
            "name": str(name).strip(),
            "is_solvent": int(code == SOLVENT),
            "charge_max": qmax,
            "charge": None,
            "charge_src": None,
            "categories": cats,
        })
    return out


def _count_pka(ws, b: "Built") -> dict:
    """{código do componente: nº de pKa} pela aba pKa (metais e sistemas fictícios ignorados)."""
    by_name = {}
    for c in b.components:
        by_name.setdefault(_norm(c["name"]), c["code"])
    out = {}
    if ws is None:
        return out
    for row in ws.iter_rows(min_row=2, values_only=True):
        name = row[0]
        if not name:
            continue
        key = _norm(name)
        if key in SKIP_SYSTEMS or METAL_RE.search(str(name)):
            continue
        code = by_name.get(ALIASES.get(key, key))
        if code is None:
            continue
        out[code] = sum(1 for x in row[2:] if _num(x) is not None)
    return out


def _assign_charges(b: "Built", aq_recs: list, n_pka: dict | None = None) -> None:
    """z_j = q_max − n_H (ver docstring do módulo)."""
    n_h: dict = dict(n_pka or {})
    for r in aq_recs:
        d = dict(r["nu"])
        if len(d) == 2 and d.get(PROTON, 0) > 0:
            (x,) = set(d) - {PROTON}
            if d[x] == 1:
                n_h[x] = max(n_h.get(x, 0), d[PROTON])
    for c in b.components:
        code, name, q = c["code"], c["name"], c["charge_max"]
        if code == SOLVENT:
            c["charge"], c["charge_src"] = 0, "solvente"
            continue
        if code == PROTON:
            c["charge"], c["charge_src"] = 1, "H+"
            if q is not None and q != 1:
                b.log.append(("alerta", f"H+ ({PROTON}): carga máxima = {q:+d}; usado +1"))
            continue
        nh = 0 if (METAL_RE.search(name) or CHARGE_PLUS_RE.search(name)) else n_h.get(code, 0)
        if q is not None:
            c["charge"], c["charge_src"] = q - nh, f"components: q_max {q:+d} − n_H {nh}"
            continue
        mr, mp = CHARGE_ROMAN_RE.search(name), CHARGE_PLUS_RE.search(name)
        if mr or mp:
            c["charge"] = ROMAN[mr.group(1)] if mr else int(mp.group(1))
            c["charge_src"] = "reserva: nome (metal)"
            b.log.append(("alerta", f"componente '{name}' ({code}): carga máxima vazia; carga deduzida "
                                    f"do nome, z = {c['charge']:+d}"))
        else:
            b.log.append(("alerta", f"componente '{name}' ({code}): carga desconhecida; preencha a "
                                    f"coluna de carga máxima da aba components (necessária p/ Debye–Hückel)"))
    for c in b.components:
        if c["charge"] is not None:
            b.log.append(("info", f"carga {c['name']} ({c['code']}): z = {c['charge']:+d} [{c['charge_src']}]"))


def _expand_categories(b: "Built") -> None:
    """'all' significa presente em todas as categorias do banco. Sem categoria: alerta."""
    real = {cat for c in b.components for cat in c["categories"] if cat != "all"}
    for c in b.components:
        cats = set(c["categories"])
        if "all" in cats:
            cats = (cats - {"all"}) | real
        if not cats:
            b.log.append(("alerta", f"componente {c['code']} ({c['name']}) sem categoria: aparece só em 'Todas'"))
        c["categories"] = sorted(cats)


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



def build(xlsx_path: Path = DEFAULT_XLSX) -> Built:
    b = Built()
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    b.components = _read_components(_sheet(wb, "components"), b)
    _expand_categories(b)

    explicit = _read_explicit_species(_sheet(wb, "species"), b)
    _hydrolysis_series_alert(explicit, b)
    aq_recs = [r for r in explicit if r["kind"] != "S0"]
    cap_recs = [r for r in explicit if r["kind"] == "S0"]
    try:
        ws_pka = _sheet(wb, "pka")
    except Exception:
        ws_pka = None
    _assign_charges(b, aq_recs, _count_pka(ws_pka, b))

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
            "INSERT INTO components(code, name, is_solvent, charge_max, charge) VALUES (?, ?, ?, ?, ?)",
            [(c["code"], c["name"], c["is_solvent"], c.get("charge_max"), c.get("charge"))
             for c in b.components],
        )
        conn.executemany(
            "INSERT INTO component_categories(comp_code, category) VALUES (?, ?)",
            [(c["code"], cat) for c in b.components for cat in c["categories"]],
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
    print(f"componentes: {len(b.components)} | espécies: {len(b.species)} | alertas/erros: {n_alert}")
    for level, msg in b.log:
        print(f"[{level}] {msg}")


if __name__ == "__main__":
    main()
