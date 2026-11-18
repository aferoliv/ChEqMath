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


def list_categories(conn: sqlite3.Connection) -> list[str]:
    rows = conn.execute("SELECT DISTINCT category FROM component_categories ORDER BY category").fetchall()
    return [r["category"] for r in rows]


def list_components(conn: sqlite3.Connection, categories: list[str] | None = None) -> list[dict]:
    """Componentes em ordem alfabética, sem repetição.

    categories=None ou vazio: todos. Caso contrário: união das categorias escolhidas
    (um componente com várias categorias aparece uma única vez).
    """
    if not categories:
        rows = conn.execute(
            "SELECT code, name, is_solvent, charge_max, charge FROM components ORDER BY name COLLATE NOCASE"
        ).fetchall()
    else:
        marks = ",".join("?" for _ in categories)
        rows = conn.execute(
            f"SELECT DISTINCT c.code, c.name, c.is_solvent, c.charge_max, c.charge FROM components c "
            f"JOIN component_categories cc ON cc.comp_code = c.code "
            f"WHERE cc.category IN ({marks}) ORDER BY c.name COLLATE NOCASE",
            categories,
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


@lru_cache(maxsize=1)
def component_charges() -> dict:
    """{código do componente: z} (só componentes com carga conhecida)."""
    conn = connect()
    try:
        rows = conn.execute("SELECT code, charge FROM components WHERE charge IS NOT NULL").fetchall()
        return {r["code"]: int(r["charge"]) for r in rows}
    finally:
        conn.close()
