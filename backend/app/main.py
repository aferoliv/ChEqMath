"""API de especiação. O front-end (frontend/) consome estes endpoints via /api/ (nginx)."""
from __future__ import annotations

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from . import engine
from .db import all_species, connect, list_components
from .engine import PROTON, SOLVENT

app = FastAPI(title="Especiação API", version="0.1.0")

VALID_KINDS = {"ABB", "complex", "S0"}
DEFAULT_KINDS = ["ABB", "complex", "S0"]


class SpeciateRequest(BaseModel):
    totals: dict[str, float]              # mol L⁻¹; H+ (0002) pode ser negativo
    kinds: list[str] = DEFAULT_KINDS
    species_codes: list[str] | None = None
    fixed_ph: float | None = None         # se informado, H+ é fixado pelo pH


class CurveRequest(BaseModel):
    totals: dict[str, float]              # ligantes/metais (sem H+)
    kinds: list[str] = DEFAULT_KINDS
    ph_min: float = 0.0
    ph_max: float = 14.0
    ph_points: int = Field(281, ge=2, le=2001)
    ref: str | None = None                # componente de referência para as frações


class MCRequest(BaseModel):
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


@app.get("/api/components")
def api_components():
    conn = connect()
    try:
        return {"items": list_components(conn)}
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
    r = _run(lambda: engine.solve(sp, totals, fixed=fixed))
    return {
        "ph": r.ph,
        "free": r.free,
        "species": r.species,
        "solids": r.solids,
        "implied_totals": r.implied_totals,
        "fractions": r.fractions,
        "residual": r.residual,
        "iterations": r.iterations,
        "active_caps": r.active_caps,
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
    data = _run(lambda: engine.ph_curve(sp, totals, ph_values, ref=ref))
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
    res = _run(lambda: engine.monte_carlo(sp, totals, fixed=fixed, sd=req.sd, n=req.n, seed=req.seed))
    res["names"] = _names(sp)
    return res
