"""Testes da correção de atividade (Debye–Hückel estendida) contra cálculos independentes."""
import math
from pathlib import Path

import pytest

from app import engine
from app.engine import PROTON, SOLVENT, ActivityModel, Species, debye_huckel_log_gamma
from app.importer import build

XLSX = Path(__file__).resolve().parents[2] / "data" / "componentes_e_especies_ChEq_Math.xlsx"
ACETATE = "0003"
CD = "0019"
A, B = engine.DH_A, engine.DH_B


@pytest.fixture(scope="module")
def built():
    return build(XLSX)


@pytest.fixture(scope="module")
def species(built):
    return [
        Species(code=s["code"], name=s["name"], kind=s["kind"], log_beta=s["log_beta"],
                nu=tuple(s["nu"]), cap_of=s["cap_of"], water_ref=bool(s.get("water_ref", 0)))
        for s in built.species
    ]


@pytest.fixture(scope="module")
def charges(built):
    return {c["code"]: c["charge"] for c in built.components if c["charge"] is not None}


def _acetic(species):
    return engine.select_species(species, {PROTON, ACETATE}, kinds=["ABB"])


def _pk(sp, nu):
    s = next(x for x in sp if x.nu_map() == nu)
    return engine.effective_log_beta(s)


# --- importador: cargas dos componentes -----------------------------------------

def test_component_charges(charges):
    assert charges[SOLVENT] == 0
    assert charges[PROTON] == 1
    assert charges[ACETATE] == -1      # q_max 0 (HAc) − 1
    assert charges[CD] == 2            # q_max +2, metal
    assert charges["0017"] == 2        # Ca(2+)
    assert charges["0018"] == -2       # carbonato: q_max 0 (H2CO3) − 2
    assert charges["0006"] == -1       # alanina: q_max +1 (H2Ala+) − 2


# --- função de Debye–Hückel -----------------------------------------------------

def test_dh_function_values():
    I = 0.01
    expected = -A * 1 * 0.1 / (1 + B * 4.5 * 0.1)
    assert debye_huckel_log_gamma(1, I, 4.5) == pytest.approx(expected, rel=1e-12)
    assert debye_huckel_log_gamma(2, I, 4.5) == pytest.approx(4 * expected, rel=1e-12)
    assert debye_huckel_log_gamma(0, I, 4.5) == 0.0
    # limite de Debye–Hückel para I → 0
    assert debye_huckel_log_gamma(1, 1e-8, 4.5) == pytest.approx(-A * 1e-4, rel=1e-3)


# --- motor --------------------------------------------------------------------------

def test_ideal_model_unchanged(species):
    sp = _acetic(species)
    r0 = engine.solve(sp, {PROTON: 0.1, ACETATE: 0.1})
    r1 = engine.solve(sp, {PROTON: 0.1, ACETATE: 0.1}, activity=ActivityModel())
    assert r1.ph == r0.ph
    assert r1.species == r0.species
    assert all(v == 0.0 for v in r1.log_gamma.values())


def test_missing_charge_raises(species):
    with pytest.raises(ValueError, match="carga desconhecida"):
        engine.solve(_acetic(species), {PROTON: 0.1, ACETATE: 0.1},
                     activity=ActivityModel(model="dh_ext", charges={PROTON: 1}))


def exact_acetic_dh(C, pKa, pKw, a=4.5):
    """HAc C mol/L com T_H+ = C, γ por DH estendida; solução independente por bissecção + ponto fixo.

    Concentrações: h + [HAc] − [OH⁻] = C;  [HAc] = C h γ² / (Ka + h γ²);  [OH⁻] = Kw / (h γ²).
    I = ½([H⁺] + [Ac⁻] + [OH⁻]) (só íons monovalentes).
    """
    Ka, Kw = 10.0 ** -pKa, 10.0 ** -pKw
    I = 0.0
    for _ in range(200):
        g = 10.0 ** debye_huckel_log_gamma(1, I, a)
        g2 = g * g

        def f(h):
            return h + C * h * g2 / (Ka + h * g2) - Kw / (h * g2) - C

        lo, hi = -14.0, 1.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(10.0 ** mid) > 0:
                hi = mid
            else:
                lo = mid
        h = 10.0 ** (0.5 * (lo + hi))
        ac = C * Ka / (Ka + h * g2)
        oh = Kw / (h * g2)
        I_new = 0.5 * (h + ac + oh)
        if abs(I_new - I) < 1e-15:
            break
        I = I_new
    return -math.log10(g * h), I


def test_acetic_acid_dh_matches_independent(species, charges):
    sp = _acetic(species)
    pKa = _pk(sp, {ACETATE: 1, PROTON: 1})
    pKw = -_pk(sp, {SOLVENT: 1, PROTON: -1})
    ph_exp, I_exp = exact_acetic_dh(0.1, pKa, pKw)
    act = ActivityModel(model="dh_ext", charges=charges)
    r = engine.solve(sp, {PROTON: 0.1, ACETATE: 0.1}, activity=act)
    assert r.ionic_strength == pytest.approx(I_exp, rel=1e-8)
    assert r.ph == pytest.approx(ph_exp, abs=1e-8)
    r0 = engine.solve(sp, {PROTON: 0.1, ACETATE: 0.1})
    # γ < 1 aumenta a dissociação em concentração ([H⁺] e [Ac⁻] maiores), mas a_H+ quase não
    # muda: em 1ª ordem a_H = √(Ka·C), independente de γ.
    assert r.free[PROTON] > r0.free[PROTON]
    assert abs(r.ph - r0.ph) < 1e-3
    assert r.residual < 1e-10


def test_fixed_ionic_strength_conditional_constant(species, charges):
    """Com I fixa e a_H+ fixa: [Ac⁻]/[HAc] = Ka / (a_H · γ_Ac)."""
    sp = _acetic(species)
    hac = next(s for s in sp if s.nu_map() == {ACETATE: 1, PROTON: 1})
    pKa = _pk(sp, {ACETATE: 1, PROTON: 1})
    I = 0.1
    act = ActivityModel(model="dh_ext", charges=charges, ionic_strength=I)
    g_ac = 10.0 ** debye_huckel_log_gamma(-1, I, 4.5)
    for pH in (3.0, 4.76, 6.0):
        r = engine.solve(sp, {ACETATE: 1e-3}, fixed={PROTON: 10.0 ** -pH}, activity=act)
        ratio = r.free[ACETATE] / r.species[hac.code]
        assert ratio == pytest.approx(10.0 ** -pKa / (10.0 ** -pH * g_ac), rel=1e-9)
        assert r.ph == pytest.approx(pH, abs=1e-12)        # pH fixado é de atividade
        assert r.ionic_strength == I


def test_ionic_strength_is_consistent(species, charges):
    sp = engine.select_species(species, {PROTON, ACETATE, CD}, kinds=["ABB", "complex"])
    act = ActivityModel(model="dh_ext", charges=charges)
    r = engine.solve(sp, {PROTON: 2e-3, ACETATE: 5e-2, CD: 1e-2}, activity=act)
    z = {s.code: sum(nu * charges[j] for j, nu in s.nu if j != SOLVENT) for s in sp}
    I = 0.5 * (sum(r.free[j] * charges[j] ** 2 for j in r.free)
               + sum(c * z[k] ** 2 for k, c in r.species.items()))
    assert r.ionic_strength == pytest.approx(I, rel=1e-8)
    for k, c in r.species.items():
        assert r.log_gamma[k] == pytest.approx(debye_huckel_log_gamma(z[k], I, 4.5), rel=1e-6, abs=1e-12)
    assert r.implied_totals[CD] == pytest.approx(1e-2, rel=1e-9)


def test_warning_above_validity_limit(species, charges):
    sp = _acetic(species)
    act = ActivityModel(model="dh_ext", charges=charges, ionic_strength=0.5)
    r = engine.solve(sp, {ACETATE: 1e-3}, fixed={PROTON: 1e-5}, activity=act)
    assert r.warnings
