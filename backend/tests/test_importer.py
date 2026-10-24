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
    assert len(built.components) == 90
    assert sum(1 for c in built.components if c["is_solvent"]) == 1
    assert len({c["code"] for c in built.components}) == len(built.components)


def test_species_ids_unique(built):
    codes = [s["code"] for s in built.species]
    assert len(codes) == len(set(codes))


def test_hydroxide_value_loaded_as_in_sheet(built):
    # Constante usada como está na planilha (coluna log K). Erro é corrigido na planilha.
    oh = next(s for s in built.species if s["name"] == "hydroxyl")
    assert oh["log_beta"] == pytest.approx(-15.75)
    assert oh["log_beta_raw"] == pytest.approx(-15.75)


def test_cadmium_series_loaded_as_in_sheet(built):
    # Coluna log K da planilha (corrigida): Cd(OH)2 cumulativo = -20,35
    rec = next(s for s in built.species if s["name"] == "Cd(OH)2(0)" and s["kind"] == "ABB")
    assert rec["log_beta"] == pytest.approx(-20.35)
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


def test_alanine_protonation_from_pka(built):
    # pKa = 2,348 e 9,867 (ordem crescente): log β1 = 9,867; log β2 = 12,215
    h1 = next(s for s in built.species if dict(s["nu"]) == {"0006": 1, PROTON: 1})
    h2 = next(s for s in built.species if dict(s["nu"]) == {"0006": 1, PROTON: 2})
    assert h1["log_beta"] == pytest.approx(9.867)
    assert h2["log_beta"] == pytest.approx(12.215)


def test_metal_hydrolysis_not_generated_from_pka(built):
    # Ca(2+) tem pKa na aba pKa, mas não deve gerar espécies protonadas
    assert not any("0017" in dict(s["nu"]) for s in built.species if s["source"].startswith("gerada"))


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
