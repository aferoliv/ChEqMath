"""Acesso de leitura ao banco de especiação (SQLite; gerado por app/importer.py).

Este módulo não lê a planilha: só consulta o banco e devolve objetos Species para o motor.
"""
from __future__ import annotations

import os
import sqlite3
from functools import lru_cache
from pathlib import Path

from .engine import Species

REPO = Path(__file__).resolve().parents[2]


def db_path() -> Path:
    return Path(os.environ.get("DB_PATH", str(REPO / "data" / "speciation.sqlite")))


def connect(path: Path | None = None) -> sqlite3.Connection:
    p = path or db_path()
    if not p.exists():
        raise FileNotFoundError(f"banco não encontrado: {p}. Rode: python -m app.importer")
    conn = sqlite3.connect(p)
    conn.row_factory = sqlite3.Row
    return conn


def list_components(conn: sqlite3.Connection) -> list[dict]:
    rows = conn.execute(
        "SELECT code, name, smiles, is_solvent FROM components ORDER BY name COLLATE NOCASE"
    ).fetchall()
    return [dict(r) for r in rows]


def load_species(conn: sqlite3.Connection) -> list[Species]:
    stoich: dict[str, list] = {}
    for r in conn.execute("SELECT species_code, comp_code, nu FROM stoich"):
        stoich.setdefault(r["species_code"], []).append((r["comp_code"], r["nu"]))
    out = []
    for r in conn.execute("SELECT code, name, kind, log_beta, cap_of, water_ref FROM species ORDER BY code"):
        out.append(
            Species(
                code=r["code"],
                name=r["name"],
                kind=r["kind"],
                log_beta=float(r["log_beta"]),
                nu=tuple(stoich.get(r["code"], [])),
                cap_of=r["cap_of"],
                water_ref=bool(r["water_ref"]),
            )
        )
    return out


@lru_cache(maxsize=1)
def all_species() -> tuple:
    conn = connect()
    try:
        return tuple(load_species(conn))
    finally:
        conn.close()
