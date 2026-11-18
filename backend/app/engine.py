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

Atividade (opcional, ver ActivityModel)
    a_i = γ_i · c_i, com γ_i pela equação de Debye–Hückel estendida:
        log10 γ_i = −A · z_i² · √I / (1 + B · å_i · √I)
        I = ½ Σ c_k z_k²  (componentes livres + espécies aquosas; sólidos e água excluídos)
    As constantes β são termodinâmicas (I = 0). Com γ fixo, a lei de ação das massas em
    concentração fica
        ln c_i = ln β_i + Σ_j ν_ij (ln [C_j] + ln γ_j) − ln γ_i
    e o Newton (em ln [C_j]) é o mesmo do caso ideal, com ln β_i corrigido. A força iônica é
    resolvida por iteração externa de ponto fixo (I → γ → solução → I), ou fixada pelo
    usuário (meio iônico).

Unidades: concentrações em mol L⁻¹; pH = −log10(a_H+). Valores em `fixed` são atividades.
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
# Debye–Hückel, água a 25 °C (Truesdell & Jones, 1974; Helgeson, 1969)
DH_A = 0.5091      # kg^½ mol^-½  (≈ L^½ mol^-½ em solução diluída)
DH_B = 0.3283      # kg^½ mol^-½ Å⁻¹
DH_ION_SIZE = 4.5  # Å, å padrão quando não há valor específico
DH_I_LIMIT = 0.1   # mol L⁻¹, limite usual de validade da DH estendida
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


@dataclass(frozen=True)
class ActivityModel:
    """Modelo de coeficientes de atividade.

    model         : 'ideal' (γ = 1) ou 'dh_ext' (Debye–Hückel estendida).
    charges       : {código do componente: z}. A carga de cada espécie é Σ ν_ij z_j.
    ion_size      : å padrão (Å).
    ion_sizes     : {código de componente ou espécie: å (Å)}, sobrepõe o padrão.
    ionic_strength: I fixa (mol L⁻¹, meio iônico). None = calculada da solução.
    """
    model: str = "ideal"
    charges: dict = field(default_factory=dict)
    ion_size: float = DH_ION_SIZE
    ion_sizes: dict = field(default_factory=dict)
    ionic_strength: Optional[float] = None
    A: float = DH_A
    B: float = DH_B

    @property
    def active(self) -> bool:
        return self.model == "dh_ext"

    def log10_gamma(self, z: np.ndarray, size: np.ndarray, I: float) -> np.ndarray:
        if not self.active or I <= 0.0:
            return np.zeros_like(z, dtype=float)
        sI = math.sqrt(I)
        return -self.A * z * z * sI / (1.0 + self.B * size * sI)


IDEAL = ActivityModel()


def debye_huckel_log_gamma(z: float, I: float, size: float = DH_ION_SIZE,
                           A: float = DH_A, B: float = DH_B) -> float:
    """log10 γ pela Debye–Hückel estendida (função escalar, útil para conferência)."""
    sI = math.sqrt(max(I, 0.0))
    return -A * z * z * sI / (1.0 + B * size * sI)


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
    converged: bool = True
    message: str = ""
    dissolved: dict = field(default_factory=dict)  # total dissolvido por componente (sem sólido)
    ionic_strength: float = 0.0        # I (mol L⁻¹), calculada ou fixada
    log_gamma: dict = field(default_factory=dict)  # log10 γ por componente e espécie aquosa
    activity_model: str = "ideal"
    I_iterations: int = 0
    warnings: list = field(default_factory=list)


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


def solve(species, totals, fixed=None, tol=1e-11, maxit=200, max_outer=30,
          activity: Optional[ActivityModel] = None, I_tol=1e-10, max_I_iter=100) -> Result:
    """Resolve o equilíbrio para os totais dados.

    totals  : {código: T_j} em mol L⁻¹ (H+ pode ser negativo).
    fixed   : {código: atividade} para componentes com a dimensão fixada (ex.: H+ por pH).
    activity: modelo de atividade (padrão: ideal, γ = 1).
    Lança ValueError (entrada inválida) ou RuntimeError (não convergiu).
    """
    activity = activity or IDEAL
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

    # --- atividade: cargas e tamanhos iônicos -------------------------------------
    lnB0, lnfix0, lnSo0 = lnB.copy(), lnfix.copy(), lnSo.copy()
    use_act = activity.active
    if use_act:
        missing = [j for j in comps if j not in activity.charges]
        if missing:
            raise ValueError(f"Debye–Hückel: carga desconhecida para o(s) componente(s) {missing}")
        if activity.ionic_strength is not None and activity.ionic_strength < 0:
            raise ValueError("força iônica fixa deve ser ≥ 0")
    zc = np.array([float(activity.charges.get(j, 0)) for j in comps])
    zq = Nq @ zc
    sz = activity.ion_sizes
    ac = np.array([float(sz.get(j, activity.ion_size)) for j in comps])
    aq_size = np.array([float(sz.get(s.code, activity.ion_size)) for s in aq])
    lnfc = np.zeros(len(comps))
    lnfq = np.zeros(len(aq))

    def apply_gamma(I):
        """Atualiza ln β, ln[C]_fixo e ln So para γ(I); devolve (lnfc, lnfq)."""
        nonlocal lnB, lnfix, lnSo
        fc = activity.log10_gamma(zc, ac, I) * LN10
        fq = activity.log10_gamma(zq, aq_size, I) * LN10
        lnB = lnB0 + Nq @ fc - fq
        lnfix = lnfix0 - fc[fidx] if fidx.size else lnfix0
        lnSo = lnSo0 - fq[tgt] if tgt.size else lnSo0
        return fc, fq

    def ionic_strength(x, c):
        return 0.5 * (float(np.sum(x * zc * zc)) + float(np.sum(c * zq * zq)))

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
        flags["ok"] = False  # devolve o último iterado; o chamador decide o que fazer
        return u, m, res, maxit

    def active_set(u, m, act):
        its = 0
        for _ in range(max_outer):
            flags["ok"] = True
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
        return u, m, act, its

    # Conjunto ativo de sólidos: começa vazio; adiciona S0 supersaturados e
    # remove os que ficam com m < 0. Com Debye–Hückel, laço externo em I.
    act: list = []
    flags = {"ok": True}
    m = np.zeros(0)
    u = np.log(np.maximum(np.abs(Tu), SCALE_FLOOR)) if Tu.size else np.zeros(0)
    I = 0.0
    I_its = 0
    if use_act and activity.ionic_strength is not None:
        I = float(activity.ionic_strength)
        lnfc, lnfq = apply_gamma(I)
        u, m, act, its = active_set(u, m, act)
        I_its = 1
    elif use_act:
        u, m, act, its = active_set(u, m, act)          # 1ª passada: γ = 1
        for I_its in range(1, max_I_iter + 1):
            lnx = lnx_full(u)
            I_new = ionic_strength(np.exp(lnx), np.exp(lnB + Nq @ lnx))
            if abs(I_new - I) <= I_tol + 1e-8 * I_new:
                I = I_new
                break
            I = I_new
            lnfc, lnfq = apply_gamma(I)
            u, m, act, its = active_set(u, m, act)
        else:
            raise RuntimeError(f"força iônica não convergiu em {max_I_iter} iterações (I = {I:.4g} mol/L)")
    else:
        u, m, act, its = active_set(u, m, act)

    F, G, c, x, _ = evaluate(u, m, act)
    res = _maxabs(scaled(F, G))
    implied = x + Nq.T @ c
    if act:
        implied = implied + Nc[np.array(act, dtype=int)].T @ m

    free = {j: float(x[cidx[j]]) for j in comps}
    species_c = {s.code: float(c[i]) for i, s in enumerate(aq)}
    solids = {caps[k].code: float(mv) for k, mv in zip(act, m)}
    impl = {j: float(implied[cidx[j]]) for j in comps}
    ph = None
    if PROTON in free:
        ph = -(math.log10(free[PROTON]) + lnfc[cidx[PROTON]] / LN10)   # pH = −log a_H+
    if not use_act:
        I = ionic_strength(x, c)  # informativa (γ = 1)
    log_gamma = {j: float(lnfc[cidx[j]] / LN10) for j in comps}
    log_gamma.update({s.code: float(lnfq[i] / LN10) for i, s in enumerate(aq)})
    warnings = []
    if use_act and I > DH_I_LIMIT:
        warnings.append(f"I = {I:.3g} mol/L > {DH_I_LIMIT} mol/L: Debye–Hückel estendida fora da faixa usual")

    # Distribuição só com espécies dissolvidas: denominador = total dissolvido de cada componente
    # (componente livre + espécies aquosas), sem o sólido precipitado.
    dissolved = {j: float(x[cidx[j]] + np.dot(Nq[:, cidx[j]], c)) for j in comps}
    fractions = {}
    for j in comps:
        if j in fixed or float(totals.get(j, 0.0)) == 0.0:
            continue
        Td = dissolved[j]
        if Td <= 0.0:
            continue
        cj = cidx[j]
        fractions[j] = {s.code: float(Nq[i, cj] * c[i] / Td)
                        for i, s in enumerate(aq) if Nq[i, cj] != 0}

    return Result(
        ph=ph,
        free=free,
        species=species_c,
        solids=solids,
        implied_totals=impl,
        fractions=fractions,
        dissolved=dissolved,
        residual=res,
        iterations=its,
        active_caps=[caps[k].code for k in act],
        converged=flags["ok"],
        message="" if flags["ok"] else f"não convergiu em {maxit} iterações (resíduo {res:.3e})",
        ionic_strength=float(I),
        log_gamma=log_gamma,
        activity_model=activity.model,
        I_iterations=I_its,
        warnings=warnings,
    )


def ph_curve(species, totals, ph_values, ref=None, activity: Optional[ActivityModel] = None):
    """Varre o pH com a_H+ fixa. `totals` são os ligantes/metais (sem H+).

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
        "free": [],  # fração do componente livre (ref) em relação ao total
        "solids": {c: [] for c in cap_codes},
        "H_total": [],
        "residual": [],
        "ionic_strength": [],
    }
    for ph in ph_values:
        r = solve(species, totals, fixed={PROTON: 10.0 ** (-float(ph))}, activity=activity)
        if not r.converged:
            raise RuntimeError(r.message)
        for c in aq_codes:
            out["species"][c].append(r.species.get(c, 0.0))
            out["fractions"][c].append(r.fractions.get(ref, {}).get(c, 0.0))
        td = r.dissolved.get(ref, 0.0)
        out["free"].append(r.free.get(ref, 0.0) / td if td > 0 else 0.0)
        for c in cap_codes:
            out["solids"][c].append(r.solids.get(c, 0.0))
        out["H_total"].append(r.implied_totals.get(PROTON, float("nan")))
        out["residual"].append(r.residual)
        out["ionic_strength"].append(r.ionic_strength)
    return out


def monte_carlo(species, totals, fixed=None, sd=0.05, n=300, seed=0,
                activity: Optional[ActivityModel] = None):
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
            r = solve(pert, totals, fixed, activity=activity)
            if not r.converged:
                raise RuntimeError(r.message)
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
