"""Testes do importador: a planilha vira componentes, espécies e estequiometria consistentes."""
from pathlib import Path

import pytest

from app.engine import PROTON, SOLVENT
from app.importer import build, write_db

XLSX = Path(__file__).resolve().parents[2] / "data" / "componentes_e_especies_ChEq_Math.xlsx"


@pytest.fixture(scope="module")
def built():
    return build(XLSX)


def test_counts(built):
    assert len(built.components) == 93
    assert sum(1 for c in built.components if c["is_solvent"]) == 1
    assert len({c["code"] for c in built.components}) == len(built.components)


def test_species_ids_unique(built):
    codes = [s["code"] for s in built.species]
    assert len(codes) == len(set(codes))


def test_hydroxide_value_loaded_as_in_sheet(built):
    # Constante usada como está na planilha (coluna log K). Erro é corrigido na planilha.
    oh = next(s for s in built.species if s["name"] == "hydroxyl")
    assert oh["log_beta"] == pytest.approx(oh["log_beta_raw"])
    assert oh["log_beta"] == pytest.approx(-15.75, abs=0.01)


def test_cadmium_series_loaded_as_in_sheet(built):
    # Coluna log K da planilha (com água incorporada): Cd(OH)2 = -23,84 (a_água = 1 → -20,35)
    rec = next(s for s in built.species if s["name"] == "Cd(OH)2(0)" and s["kind"] == "ABB")
    assert rec["log_beta"] == pytest.approx(-23.84, abs=0.01)
    assert any(lvl == "alerta" and "Cd(OH)2(0)" in msg for lvl, msg in built.log)


def test_comma_decimal_parsed(built):
    s0 = next(s for s in built.species if s["kind"] == "S0")
    assert s0["log_beta"] == pytest.approx(-6.7)


def test_s0_links_to_neutral_aqueous_species(built):
    s0 = next(s for s in built.species if s["kind"] == "S0")
    target = next(s for s in built.species if s["code"] == s0["cap_of"])
    assert target["kind"] != "S0"
    assert target["name"] == "Cd(OH)2(0)"
    assert dict(target["nu"]) == dict(s0["nu"])


def test_acetic_acid_not_duplicated(built):
    hits = [s for s in built.species if dict(s["nu"]) == {"0003": 1, PROTON: 1}]
    assert len(hits) == 1
    assert hits[0]["source"] == "planilha:Species"


def test_pka_sheet_not_used(built):
    # A aba pKa é apenas de apoio: nenhuma espécie é gerada a partir dela.
    assert all(s["source"] == "planilha:Species" for s in built.species)


def test_all_stoich_codes_exist(built):
    comp_codes = {c["code"] for c in built.components}
    for s in built.species:
        for j, _ in s["nu"]:
            assert j in comp_codes, f"{s['name']} usa componente desconhecido {j}"


def test_water_has_no_species_row_of_its_own(built):
    assert all(s["code"] != SOLVENT for s in built.species)


def test_write_db_roundtrip(built, tmp_path):
    db = tmp_path / "t.sqlite"
    write_db(built, db)
    import sqlite3

    conn = sqlite3.connect(db)
    n_sp = conn.execute("SELECT COUNT(*) FROM species").fetchone()[0]
    n_st = conn.execute("SELECT COUNT(*) FROM stoich").fetchone()[0]
    conn.close()
    assert n_sp == len(built.species)
    assert n_st == sum(len(s["nu"]) for s in built.species)


# --- carga dos componentes: coluna B = carga máxima do sistema ----------------------

@pytest.mark.parametrize("txt, z", [
    ("2+", 2), ("+", 1), ("-", -1), ("3-", -3), ("+2", 2), ("-1", -1), ("0", 0),
    ("−", -1), (" 2 + ", 2), (2.0, 2), (-1, -1), (None, None), ("", None),
])
def test_parse_charge(txt, z):
    from app.importer import _parse_charge
    assert _parse_charge(txt) == z


@pytest.mark.parametrize("bad", ["2+-", "+-", "abc", "1.5", 1.5])
def test_parse_charge_invalid(bad):
    from app.importer import _parse_charge
    with pytest.raises(ValueError):
        _parse_charge(bad)


def test_charge_from_charge_column(tmp_path):
    """Coluna de carga máxima (achada pelo cabeçalho) = carga da espécie totalmente protonada;
    z = q_max − n_H."""
    import openpyxl
    from app.importer import _component_columns
    wb = openpyxl.load_workbook(XLSX)
    ws = next(w for w in wb.worksheets if w.title.lower().startswith("components"))
    col = _component_columns([c.value for c in ws[1]])
    assert col["charge"] is not None
    qcol = {"0019": "2+", "0038": "+", "0003": "0", "0030": "+", "0044": "-", "0029": "2+",
            "0002": "+", "0001": "0"}
    for row in ws.iter_rows(min_row=2):
        code = row[col["code"]].value
        if code is not None and str(code).strip() in qcol:
            row[col["charge"]].value = qcol[str(code).strip()]
    path = tmp_path / "plan.xlsx"
    wb.save(path)
    b = build(path)
    z = {c["code"]: c["charge"] for c in b.components}
    assert z["0019"] == 2      # Cd(II): metal → z = q_max
    assert z["0038"] == -1     # glicina: +1 − 2 pKa
    assert z["0003"] == -1     # acetato: 0 − 1
    assert z["0030"] == 0      # etanolamina: +1 − 1
    assert z["0044"] == -1     # fluoreto: sem pKa → z = q_max
    assert z["0029"] == -4     # EDTA: +2 − 6
    assert z[PROTON] == 1 and z[SOLVENT] == 0


def test_component_columns_by_header():
    from app.importer import _component_columns
    new = ["Components", "Maximum Electrical Charge", "SMILES", "code", "category 1", "category 2", "category 3"]
    assert _component_columns(new) == {"name": 0, "charge": 1, "code": 3, "cats": [4, 5, 6]}
    old = ["Components", "SMILES", "code", "category 1", "category 2", "category 3"]
    assert _component_columns(old) == {"name": 0, "charge": None, "code": 2, "cats": [3, 4, 5]}
