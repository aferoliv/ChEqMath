"""Testes do motor de especiação contra soluções independentes."""
import math
import random
from pathlib import Path

import pytest

from app import engine
from app.engine import PROTON, SOLVENT, Species
from app.importer import build

XLSX = Path(__file__).resolve().parents[2] / "data" / "componentes_e_especies_ChEq_Math.xlsx"
ACETATE = "0003"
ALANINE = "0006"
CD = "0019"


@pytest.fixture(scope="module")
def species():
    return _to_species(build(XLSX))


def _to_species(b):
    return [
        Species(code=s["code"], name=s["name"], kind=s["kind"], log_beta=s["log_beta"],
                nu=tuple(s["nu"]), cap_of=s["cap_of"], water_ref=bool(s.get("water_ref", 0)))
        for s in b.species
    ]


def _find(species, nu_dict):
    for s in species:
        if s.nu_map() == nu_dict:
            return s
    raise KeyError(nu_dict)


def exact_ph_acetic(C, pka=4.76, pkw=14.0):
    """pH exato de HAc C mol/L com T_H+ = C (sem titulante): h + HAc − OH = C."""
    Ka, Kw = 10.0 ** -pka, 10.0 ** -pkw

    def g(h):
        return h + C * h / (Ka + h) - Kw / h - C

    lo, hi = -14.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(10.0 ** mid) > 0:
            hi = mid
        else:
            lo = mid
    return -0.5 * (lo + hi)


def bjerrum_fractions(pkas_asc, pH):
    """Frações de H_kX (k = 0..n) pela equação de Bjerrum, independentes do motor."""
    n = len(pkas_asc)
    h = 10.0 ** -pH
    betas = [1.0] + [10.0 ** sum(pkas_asc[n - k:]) for k in range(1, n + 1)]
    terms = [betas[k] * h ** k for k in range(n + 1)]
    s = sum(terms)
    return [t / s for t in terms]


def test_acetic_acid_pH_matches_exact_root(species):
    sp = engine.select_species(species, {PROTON, ACETATE}, kinds=["ABB"])
    pka = _find(sp, {ACETATE: 1, PROTON: 1}).log_beta          # valor da planilha
    pkw = -engine.effective_log_beta(_find(sp, {SOLVENT: 1, PROTON: -1}))
    r = engine.solve(sp, {PROTON: 0.1, ACETATE: 0.1})
    assert r.ph == pytest.approx(exact_ph_acetic(0.1, pka, pkw), abs=1e-8)
    assert r.residual < 1e-10


def test_monoprotic_fractions_are_bjerrum(species):
    sp = engine.select_species(species, {PROTON, ACETATE}, kinds=["ABB"])
    hac = _find(sp, {ACETATE: 1, PROTON: 1})
    pka = hac.log_beta                                           # valor da planilha
    for pH in (1.0, 3.0, 4.76, 6.0, 9.0, 13.0):
        r = engine.solve(sp, {ACETATE: 1e-3}, fixed={PROTON: 10.0 ** -pH})
        frac_hac = r.fractions[ACETATE][hac.code]
        frac_ac = r.free[ACETATE] / 1e-3
        expected = bjerrum_fractions([pka], pH)
        assert frac_ac == pytest.approx(expected[0], abs=1e-9)
        assert frac_hac == pytest.approx(expected[1], abs=1e-9)


@pytest.mark.skip(reason="alanina sem espécies na aba Species; aba pKa não é usada")
def test_diprotic_alanine_fractions_are_bjerrum(species):
    sp = engine.select_species(species, {PROTON, ALANINE}, kinds=["ABB"])
    h1 = _find(sp, {ALANINE: 1, PROTON: 1})
    h2 = _find(sp, {ALANINE: 1, PROTON: 2})
    pkas = [2.348, 9.867]  # valores da planilha (aba pKa)
    for pH in (1.0, 2.35, 6.0, 9.87, 12.0):
        r = engine.solve(sp, {ALANINE: 0.01}, fixed={PROTON: 10.0 ** -pH})
        exp = bjerrum_fractions(pkas, pH)
        assert r.free[ALANINE] / 0.01 == pytest.approx(exp[0], abs=1e-9)
        assert r.fractions[ALANINE][h1.code] == pytest.approx(exp[1], abs=1e-9)
        assert r.fractions[ALANINE][h2.code] == pytest.approx(exp[2], abs=1e-9)


def test_fractions_sum_to_one_with_solid(species):
    sp = engine.select_species(species, {PROTON, CD}, kinds=["ABB", "S0"])
    r = engine.solve(sp, {CD: 1e-2}, fixed={PROTON: 1e-12})
    # distribuição só com espécies dissolvidas: soma = 1 sobre o total dissolvido
    total = r.free[CD] / r.dissolved[CD] + sum(r.fractions[CD].values())
    assert total == pytest.approx(1.0, abs=1e-9)
    assert r.solids and all(code not in r.fractions[CD] for code in r.solids)


def test_cd_hydroxide_solubility_cap(species):
    sp = engine.select_species(species, {PROTON, CD}, kinds=["ABB", "S0"])
    cap = next(s for s in sp if s.kind == "S0")
    aq = next(s for s in sp if s.code == cap.cap_of)
    So = 10.0 ** cap.log_beta

    # Varredura de pH com Cd total alto: em algum pH o Cd(OH)2 aquoso atinge So
    # e o excesso precipita. Em pH 12 o Cd(OH)4^2- domina e o limite não ativa.
    binding = []
    for pH in [6.0 + 0.25 * i for i in range(25)]:
        r = engine.solve(sp, {CD: 1e-2}, fixed={PROTON: 10.0 ** -pH})
        if cap.code in r.active_caps:
            binding.append((pH, r))
    assert binding, "o limite de solubilidade deve ativar em alguma faixa de pH"
    for pH, r in binding:
        assert r.species[aq.code] == pytest.approx(So, rel=1e-8)
        assert r.solids[cap.code] > 0
        assert r.implied_totals[CD] == pytest.approx(1e-2, rel=1e-8)

    # Total pequeno em pH ácido: sem sólido, espécie abaixo de So
    r2 = engine.solve(sp, {CD: 1e-7}, fixed={PROTON: 1e-4})
    assert r2.active_caps == []
    assert r2.species[aq.code] < So


def test_random_mass_balance_closes(species):
    sp = engine.select_species(species, {PROTON, ALANINE, CD}, kinds=["ABB", "complex", "S0"])
    rng = random.Random(7)
    ok = 0
    for _ in range(12):
        totals = {
            PROTON: rng.uniform(-2e-3, 2e-3),
            ALANINE: rng.uniform(1e-4, 1e-2),
            CD: rng.uniform(1e-5, 1e-3),
        }
        try:
            r = engine.solve(sp, totals)
        except (RuntimeError, ValueError):
            continue
        ok += 1
        for j, T in totals.items():
            assert r.implied_totals[j] == pytest.approx(T, abs=1e-9)
    assert ok >= 8, "solver deve convergir na maioria dos casos aleatórios"


def test_missing_component_total_raises(species):
    sp = engine.select_species(species, {PROTON, ACETATE}, kinds=["ABB"])
    with pytest.raises(ValueError):
        engine.solve(sp, {PROTON: 0.1})  # acetato sem total


def test_monte_carlo_zero_sd_is_deterministic(species):
    sp = engine.select_species(species, {PROTON, ACETATE}, kinds=["ABB"])
    det = engine.solve(sp, {PROTON: 0.1, ACETATE: 0.1}).ph
    mc = engine.monte_carlo(sp, {PROTON: 0.1, ACETATE: 0.1}, sd=0.0, n=5)
    assert mc["pH"]["mean"] == pytest.approx(det, abs=1e-9)
    assert mc["pH"]["std"] == pytest.approx(0.0, abs=1e-9)


def test_monte_carlo_spread_grows_with_sd(species):
    sp = engine.select_species(species, {PROTON, ACETATE}, kinds=["ABB"])
    small = engine.monte_carlo(sp, {PROTON: 0.1, ACETATE: 0.1}, sd=0.01, n=60, seed=1)
    large = engine.monte_carlo(sp, {PROTON: 0.1, ACETATE: 0.1}, sd=0.2, n=60, seed=1)
    assert large["pH"]["std"] > small["pH"]["std"] > 0


def test_water_convention_hydroxide_gives_kw(species):
    # Constante da planilha (≈ −15,75) referida a água 55,6 mol/L → K_w ≈ 10^-14 com a_água = 1
    from app.engine import effective_log_beta, WATER_CONC
    oh = next(s for s in species if s.nu_map() == {SOLVENT: 1, PROTON: -1})
    assert effective_log_beta(oh) == pytest.approx(oh.log_beta + math.log10(WATER_CONC), abs=1e-12)
    assert effective_log_beta(oh) == pytest.approx(-14.0, abs=0.01)
    kw = 10.0 ** effective_log_beta(oh)
    sp = engine.select_species(species, {PROTON}, kinds=["ABB"])
    r = engine.solve(sp, {PROTON: 0.0}, fixed={PROTON: 1e-8})
    oh_conc = r.species[oh.code]
    assert 1e-8 * oh_conc == pytest.approx(kw, rel=1e-6)


def test_water_convention_leaves_non_water_species_unchanged(species):
    from app.engine import effective_log_beta
    ac = _find(species, {ACETATE: 1, PROTON: 1})
    assert effective_log_beta(ac) == pytest.approx(ac.log_beta, abs=1e-12)
