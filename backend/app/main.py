"""API de especiação. O front-end (frontend/) consome estes endpoints via /api/ (nginx)."""
from __future__ import annotations

from typing import Literal

from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel, Field

from . import engine
from .db import all_species, component_charges, connect, list_categories, list_components
from .engine import DH_ION_SIZE, PROTON, SOLVENT, ActivityModel

APP_VERSION = "0.3.1"   # mesma versão do pacote zip e de FRONT_VERSION em frontend/app.js
app = FastAPI(title="Especiação API", version=APP_VERSION)

VALID_KINDS = {"ABB", "complex", "S0"}
DEFAULT_KINDS = ["ABB", "complex", "S0"]


class ActivityOptions(BaseModel):
    """Opções de atividade, comuns a todos os cálculos."""
    activity: Literal["ideal", "dh_ext"] = "ideal"   # dh_ext = Debye–Hückel estendida
    ion_size: float = Field(DH_ION_SIZE, gt=0, le=20)  # å padrão (Å)
    ion_sizes: dict[str, float] | None = None          # å por código de componente/espécie
    ionic_strength: float | None = Field(None, ge=0)   # I fixa (meio iônico); None = calculada


class SpeciateRequest(ActivityOptions):
    totals: dict[str, float]              # mol L⁻¹; H+ (0002) pode ser negativo
    kinds: list[str] = DEFAULT_KINDS
    species_codes: list[str] | None = None
    fixed_ph: float | None = None         # se informado, H+ é fixado pelo pH


class CurveRequest(ActivityOptions):
    totals: dict[str, float]              # ligantes/metais (sem H+)
    kinds: list[str] = DEFAULT_KINDS
    ph_min: float = 0.0
    ph_max: float = 14.0
    ph_points: int = Field(281, ge=2, le=2001)
    ref: str | None = None                # componente de referência para as frações


class MCRequest(ActivityOptions):
    totals: dict[str, float]
    kinds: list[str] = DEFAULT_KINDS
    fixed_ph: float | None = None
    sd: float = Field(0.05, ge=0.0)       # desvio-padrão de log β (unidades de log10)
    n: int = Field(300, ge=10, le=5000)
    seed: int = 0


def _check_kinds(kinds: list[str]) -> list[str]:
    bad = set(kinds) - VALID_KINDS
    if bad:
        raise HTTPException(422, f"tipos de equilíbrio inválidos: {sorted(bad)}")
    return kinds


def _totals(raw: dict[str, float]) -> dict[str, float]:
    """Descarta totais nulos e a água (solvente)."""
    return {k: float(v) for k, v in raw.items() if k != SOLVENT and float(v) != 0.0}


def _select(totals: dict[str, float], kinds: list[str], codes: list[str] | None = None):
    allowed = set(totals) | {PROTON}
    return engine.select_species(all_species(), allowed, kinds=kinds, codes=codes)


def _activity(req: ActivityOptions) -> ActivityModel:
    if req.activity == "ideal":
        return ActivityModel()
    return ActivityModel(
        model="dh_ext",
        charges=component_charges(),
        ion_size=req.ion_size,
        ion_sizes=dict(req.ion_sizes or {}),
        ionic_strength=req.ionic_strength,
    )


def _names(species) -> dict[str, str]:
    return {s.code: s.name for s in species}


def _run(fn):
    try:
        return fn()
    except (ValueError, RuntimeError) as exc:
        raise HTTPException(422, str(exc))


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.get("/api/version")
def api_version():
    return {"version": APP_VERSION}


@app.get("/api/categories")
def api_categories():
    conn = connect()
    try:
        return {"items": list_categories(conn)}
    finally:
        conn.close()


@app.get("/api/components")
def api_components(categories: list[str] = Query(default=[])):
    """Componentes filtrados pela união das categorias (sem parâmetro: todos).
    O valor '__all__' equivale a sem filtro."""
    cats = [c for c in categories if c and c != "__all__"]
    conn = connect()
    try:
        return {"items": list_components(conn, cats or None)}
    finally:
        conn.close()


@app.get("/api/species")
def api_species(kinds: str = ",".join(DEFAULT_KINDS)):
    ks = _check_kinds(kinds.split(","))
    items = [
        {"code": s.code, "name": s.name, "kind": s.kind, "log_beta": s.log_beta,
         "nu": [list(p) for p in s.nu], "cap_of": s.cap_of}
        for s in all_species() if s.kind in ks
    ]
    return {"items": items}


@app.post("/api/speciate")
def api_speciate(req: SpeciateRequest):
    _check_kinds(req.kinds)
    totals = _totals(req.totals)
    if not totals and req.fixed_ph is None:
        raise HTTPException(422, "informe ao menos um total ou um pH fixo")
    sp = _select(totals, req.kinds, req.species_codes)
    fixed = {PROTON: 10.0 ** (-req.fixed_ph)} if req.fixed_ph is not None else None
    act = _activity(req)
    r = _run(lambda: engine.solve(sp, totals, fixed=fixed, activity=act))
    return {
        "ph": r.ph,
        "free": r.free,
        "species": r.species,
        "solids": r.solids,
        "implied_totals": r.implied_totals,
        "fractions": r.fractions,
        "residual": r.residual,
        "iterations": r.iterations,
        "converged": r.converged,
        "message": r.message,
        "active_caps": r.active_caps,
        "activity_model": r.activity_model,
        "ionic_strength": r.ionic_strength,
        "log_gamma": r.log_gamma,
        "I_iterations": r.I_iterations,
        "warnings": r.warnings,
        "names": _names(sp),
    }


@app.post("/api/curve")
def api_curve(req: CurveRequest):
    _check_kinds(req.kinds)
    totals = {k: v for k, v in _totals(req.totals).items() if k != PROTON}
    if not totals:
        raise HTTPException(422, "informe ao menos um total de ligante ou metal")
    ref = req.ref or next(iter(totals))
    if ref not in totals:
        raise HTTPException(422, f"referência {ref} não está entre os totais informados")
    sp = _select(totals, req.kinds)
    ph_values = [req.ph_min + i * (req.ph_max - req.ph_min) / (req.ph_points - 1)
                 for i in range(req.ph_points)]
    act = _activity(req)
    data = _run(lambda: engine.ph_curve(sp, totals, ph_values, ref=ref, activity=act))
    data["activity_model"] = act.model
    data["names"] = _names(sp)
    return data


@app.post("/api/mc")
def api_mc(req: MCRequest):
    _check_kinds(req.kinds)
    totals = _totals(req.totals)
    if not totals and req.fixed_ph is None:
        raise HTTPException(422, "informe ao menos um total ou um pH fixo")
    sp = _select(totals, req.kinds)
    fixed = {PROTON: 10.0 ** (-req.fixed_ph)} if req.fixed_ph is not None else None
    act = _activity(req)
    res = _run(lambda: engine.monte_carlo(sp, totals, fixed=fixed, sd=req.sd, n=req.n, seed=req.seed,
                                          activity=act))
    res["activity_model"] = act.model
    res["names"] = _names(sp)
    return res
