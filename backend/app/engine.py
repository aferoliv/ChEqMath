"""Motor de especiação (modelo componentes–espécies, Teoria.docx).

Modelo
    a_i = 10**log_beta_i · Π_j a_j**nu_ij          (I = 0: atividade ≈ concentração)
    T_j = [C_j] + Σ_i nu_ij · c_i + Σ_s nu_sj · m_s   (balanço de massa do componente j)

Sólidos (espécies S0)
    Limite de solubilidade da espécie aquosa alvo: c_alvo = So, com So = 10**log_beta_S0.
    Se a espécie aquosa supera So, o excesso precipita (m_s ≥ 0).

Convenções
    - nu_ij: unidades do componente j na espécie i (negativo quando a espécie "consome" j).
    - Água (0001) é o solvente e fica fora do sistema de incógnitas. Suas constantes
      são convertidas por f(água) = 1/55,6 (ver activity_coefficient); demais f = 1.
    - H+ (0002) é o componente de referência do pH. Seu total T_H+ pode ser negativo
      (adição de base) e, no modo de pH fixo, é calculado a partir da solução.

Unidades: concentrações em mol L⁻¹; pH = −log10(a_H+).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np

LN10 = math.log(10.0)
SOLVENT = "0001"  # água
PROTON = "0002"   # H+: referência do pH
SCALE_FLOOR = 1e-6  # piso para normalizar resíduos de balanço (mol L⁻¹)
WATER_CONC = 55.6  # mol L⁻¹, concentração molar da água usada nas constantes da planilha


def activity_coefficient(code: str) -> float:
    """Coeficiente de atividade f_j (provisório).

    f(água) = 1/55,6 e f = 1 para as demais espécies. As constantes da planilha
    referem-se à água em concentração molar (55,6 mol L⁻¹); f(água) converte-as
    para a convenção a_água = 1. Funções de atividade (Debye–Hückel etc.) entrarão aqui.
    """
    return 1.0 / WATER_CONC if code == SOLVENT else 1.0


def effective_log_beta(s) -> float:
    """log β convertido para a convenção a_água = 1.

    Só para espécies com water_ref=True (constante já referida a 55,6 mol L⁻¹):
        log β_ef = log β − ν_água · log10 f(água)
    Hidróxido: ν_água = +1 → log β_ef = −15,75 + 1,745 = −14,00 (K_w).
    Demais espécies (inclusive as de Cd, em convenção a_água = 1) não são alteradas.
    """
    if not s.water_ref:
        return s.log_beta
    nu_w = dict(s.nu).get(SOLVENT, 0)
    return s.log_beta - nu_w * math.log10(activity_coefficient(SOLVENT))


@dataclass(frozen=True)
class Species:
    code: str
    name: str
    kind: str                     # 'ABB' | 'complex' | 'S0'
    log_beta: float               # a_i = 10**log_beta · Π a_j**nu_ij
    nu: tuple                     # ((codigo_componente, nu), ...)
    cap_of: Optional[str] = None  # S0: código da espécie aquosa limitada pela solubilidade
    water_ref: bool = False       # constante referida a água 55,6 mol/L (converter para a_água = 1)

    def nu_map(self) -> dict:
        return dict(self.nu)

    def components(self) -> set:
        return {c for c, _ in self.nu if c != SOLVENT}


@dataclass
class Result:
    ph: Optional[float]
    free: dict            # [C_j] mol L⁻¹ (componentes fixos incluídos)
    species: dict         # c_i mol L⁻¹ (espécies aquosas)
    solids: dict          # m_s mol L⁻¹ (sólido precipitado), por código S0
    implied_totals: dict  # T_j recalculado a partir da solução
    fractions: dict       # comp -> {código -> fração de T_j} (espécies e sólidos)
    residual: float       # máx |balanço|/escala e |saturação| no ponto final
    iterations: int
    active_caps: list = field(default_factory=list)


def select_species(species, allowed, kinds=None, codes=None):
    """Mantém as espécies cujos componentes estão em `allowed` (a água é sempre permitida).

    Um S0 só é mantido se a espécie aquosa que ele limita também estiver selecionada.
    """
    allowed = set(allowed) | {SOLVENT}
    kinds = set(kinds) if kinds is not None else None
    codes = set(codes) if codes is not None else None
    out = []
    for s in species:
        if kinds is not None and s.kind not in kinds:
            continue
        if codes is not None and s.code not in codes:
            continue
        if s.components() <= allowed:
            out.append(s)
    aq_codes = {s.code for s in out if s.kind != "S0"}
    return [s for s in out if s.kind != "S0" or s.cap_of in aq_codes]


def _maxabs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a))) if a.size else 0.0


def solve(species, totals, fixed=None, tol=1e-11, maxit=200, max_outer=30) -> Result:
    """Resolve o equilíbrio para os totais dados.

    totals : {código: T_j} em mol L⁻¹ (H+ pode ser negativo).
    fixed  : {código: atividade} para componentes com a dimensão fixada (ex.: H+ por pH).
    Lança ValueError (entrada inválida) ou RuntimeError (não convergiu).
    """
    fixed = {j: float(v) for j, v in (fixed or {}).items()}
    if any(v <= 0 for v in fixed.values()):
        raise ValueError("valores fixos (atividades) devem ser > 0")

    aq = [s for s in species if s.kind != "S0"]
    caps = [s for s in species if s.kind == "S0"]
    aq_idx = {s.code: i for i, s in enumerate(aq)}
    for c in caps:
        if c.cap_of not in aq_idx:
            raise ValueError(f"S0 {c.code} aponta para espécie aquosa ausente ({c.cap_of})")

    comps = sorted({j for s in species for j in s.components()} | set(totals) | set(fixed))
    comps = [j for j in comps if j != SOLVENT]
    cidx = {j: k for k, j in enumerate(comps)}
    unk = [j for j in comps if j not in fixed]
    for j in unk:
        if j not in totals and j != PROTON:
            raise ValueError(f"componente {j} sem total informado")
    Tu = np.array([float(totals.get(j, 0.0)) for j in unk], dtype=float)
    scale = np.maximum(np.abs(Tu), SCALE_FLOOR)
    uidx = np.array([cidx[j] for j in unk], dtype=int)
    fidx = np.array([cidx[j] for j in fixed], dtype=int)
    lnfix = np.log(np.array([fixed[j] for j in fixed], dtype=float))

    Nq = np.zeros((len(aq), len(comps)))
    for i, s in enumerate(aq):
        for j, nu in s.nu:
            if j in cidx:
                Nq[i, cidx[j]] = nu
    lnB = np.array([effective_log_beta(s) * LN10 for s in aq], dtype=float)

    Nc = np.zeros((len(caps), len(comps)))
    for k, s in enumerate(caps):
        for j, nu in s.nu:
            if j in cidx:
                Nc[k, cidx[j]] = nu
    tgt = np.array([aq_idx[c.cap_of] for c in caps], dtype=int)
    lnSo = np.array([c.log_beta * LN10 for c in caps], dtype=float)

    def lnx_full(u):
        lnx = np.zeros(len(comps))
        lnx[uidx] = u
        lnx[fidx] = lnfix
        return lnx

    def evaluate(u, m, act):
        lnx = lnx_full(u)
        lnc = lnB + Nq @ lnx
        c = np.exp(lnc)
        x = np.exp(lnx)
        F = x[uidx] + (Nq.T @ c)[uidx] - Tu
        if act:
            A = np.array(act, dtype=int)
            F = F + (Nc[A].T @ m)[uidx]
            G = lnc[tgt[A]] - lnSo[A]
        else:
            G = np.zeros(0)
        return F, G, c, x, lnc

    def scaled(F, G):
        parts = []
        if F.size:
            parts.append(F / scale)
        if G.size:
            parts.append(G)
        return np.concatenate(parts) if parts else np.zeros(0)

    def jacobian(c, x, act):
        nu_, na = len(unk), len(act)
        J = np.zeros((nu_ + na, nu_ + na))
        if nu_:
            M = Nq.T @ (Nq * c[:, None])
            J[:nu_, :nu_] = M[np.ix_(uidx, uidx)] + np.diag(x[uidx])
        if na:
            A = np.array(act, dtype=int)
            J[:nu_, nu_:] = Nc[A][:, uidx].T
            J[nu_:, :nu_] = Nq[tgt[A]][:, uidx]
        return J

    def newton(u, m, act):
        u, m = u.copy(), m.copy()
        n_var = len(unk) + len(act)
        if n_var == 0:
            return u, m, 0.0, 0
        F, G, c, x, _ = evaluate(u, m, act)
        res = _maxabs(scaled(F, G))
        for it in range(maxit):
            if res < tol:
                return u, m, res, it
            J = jacobian(c, x, act)
            dw = -np.linalg.solve(J, np.concatenate([F, G]))
            du, dm = dw[: len(unk)], dw[len(unk):]
            lam = 1.0
            if du.size:
                big = float(np.max(np.abs(du)))
                if big > 2.0:
                    lam = 2.0 / big
            while True:
                un, mn = u + lam * du, m + lam * dm
                Fn, Gn, cn, xn, _ = evaluate(un, mn, act)
                resn = _maxabs(scaled(Fn, Gn))
                if resn < res or lam < 1e-8:
                    break
                lam *= 0.5
            u, m, F, G, c, x, res = un, mn, Fn, Gn, cn, xn, resn
        raise RuntimeError(f"não convergiu em {maxit} iterações (resíduo {res:.2e})")

    # Conjunto ativo de sólidos: começa vazio; adiciona S0 supersaturados e
    # remove os que ficam com m < 0.
    act: list = []
    m = np.zeros(0)
    u = np.log(np.maximum(np.abs(Tu), SCALE_FLOOR)) if Tu.size else np.zeros(0)
    its = 0
    for _ in range(max_outer):
        u, m, _, its = newton(u, m, act)
        new_act, new_m = [], []
        for k, mv in zip(act, m):
            if mv >= -1e-15:
                new_act.append(k)
                new_m.append(mv)
        changed = len(new_act) != len(act)
        act, m = new_act, np.array(new_m, dtype=float)
        lnc_now = lnB + Nq @ lnx_full(u)
        for k in range(len(caps)):
            if k not in act and lnc_now[tgt[k]] > lnSo[k] + 1e-9:
                act.append(k)
                m = np.append(m, 0.0)
                changed = True
        if not changed:
            break
    else:
        raise RuntimeError("conjunto de sólidos ativos não estabilizou")

    F, G, c, x, _ = evaluate(u, m, act)
    res = _maxabs(scaled(F, G))
    implied = x + Nq.T @ c
    if act:
        implied = implied + Nc[np.array(act, dtype=int)].T @ m

    free = {j: float(x[cidx[j]]) for j in comps}
    species_c = {s.code: float(c[i]) for i, s in enumerate(aq)}
    solids = {caps[k].code: float(mv) for k, mv in zip(act, m)}
    impl = {j: float(implied[cidx[j]]) for j in comps}
    ph = -math.log10(free[PROTON]) if PROTON in free else None

    fractions = {}
    for j in comps:
        if j in fixed or float(totals.get(j, 0.0)) == 0.0:
            continue
        Tj = float(totals[j])
        cj = cidx[j]
        fr = {}
        for i, s in enumerate(aq):
            if Nq[i, cj] != 0:
                fr[s.code] = float(Nq[i, cj] * c[i] / Tj)
        for k, mv in zip(act, m):
            if Nc[k, cj] != 0:
                fr[caps[k].code] = float(Nc[k, cj] * mv / Tj)
        fractions[j] = fr

    return Result(
        ph=ph,
        free=free,
        species=species_c,
        solids=solids,
        implied_totals=impl,
        fractions=fractions,
        residual=res,
        iterations=its,
        active_caps=[caps[k].code for k in act],
    )


def ph_curve(species, totals, ph_values, ref=None):
    """Varre o pH com H+ fixo. `totals` são os ligantes/metais (sem H+).

    Retorna concentrações, frações relativas ao componente `ref` e o total de H+
    implícito (prótons adicionados, mol L⁻¹, negativo para base) em cada pH.
    """
    if not totals:
        raise ValueError("informe ao menos um total")
    ref = ref or next(iter(totals))
    aq_codes = [s.code for s in species if s.kind != "S0"]
    cap_codes = [s.code for s in species if s.kind == "S0"]
    out = {
        "ph": [float(p) for p in ph_values],
        "ref": ref,
        "species": {c: [] for c in aq_codes},
        "fractions": {c: [] for c in aq_codes},
        "solids": {c: [] for c in cap_codes},
        "H_total": [],
        "residual": [],
    }
    for ph in ph_values:
        r = solve(species, totals, fixed={PROTON: 10.0 ** (-float(ph))})
        for c in aq_codes:
            out["species"][c].append(r.species.get(c, 0.0))
            out["fractions"][c].append(r.fractions.get(ref, {}).get(c, 0.0))
        for c in cap_codes:
            out["solids"][c].append(r.solids.get(c, 0.0))
        out["H_total"].append(r.implied_totals.get(PROTON, float("nan")))
        out["residual"].append(r.residual)
    return out


def monte_carlo(species, totals, fixed=None, sd=0.05, n=300, seed=0):
    """Propagação de incerteza: log β ~ N(log β, sd) em cada espécie (sd em unidades de log10).

    O sd não vem da planilha; é um parâmetro do usuário (padrão 0,05).
    """
    rng = np.random.default_rng(seed)
    ph_vals, sp_vals, fails = [], {}, 0
    for _ in range(n):
        pert = [
            replace(s, log_beta=s.log_beta + float(rng.normal(0.0, sd)))
            for s in species
        ]
        try:
            r = solve(pert, totals, fixed)
        except (RuntimeError, ValueError, FloatingPointError, OverflowError):
            fails += 1
            continue
        if r.ph is not None:
            ph_vals.append(r.ph)
        for code, val in r.species.items():
            sp_vals.setdefault(code, []).append(val)

    def stats(v):
        arr = np.asarray(v, dtype=float)
        return {
            "mean": float(arr.mean()),
            "std": float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
            "n": int(arr.size),
        }

    return {
        "n_ok": n - fails,
        "n_fail": fails,
        "sd_logbeta": sd,
        "pH": stats(ph_vals) if ph_vals else None,
        "species": {c: stats(v) for c, v in sp_vals.items()},
    }
