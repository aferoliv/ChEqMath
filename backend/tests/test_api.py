"""Testes da API (FastAPI TestClient). Usa o banco gerado pelo importador."""
import os
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
DB = REPO / "data" / "speciation.sqlite"


@pytest.fixture(scope="module")
def client():
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient

    if not DB.exists():
        from app.importer import build, write_db
        write_db(build(REPO / "data" / "componentes_e_especies_ChEq_Math.xlsx"), DB)
    os.environ.setdefault("DB_PATH", str(DB))
    from app.main import app
    return TestClient(app)


def test_components(client):
    r = client.get("/api/components")
    assert r.status_code == 200
    assert len(r.json()["items"]) == 90


def test_speciate_acetic_acid(client):
    r = client.post("/api/speciate", json={"totals": {"0002": 0.1, "0003": 0.1}, "kinds": ["ABB"]})
    assert r.status_code == 200
    assert r.json()["ph"] == pytest.approx(2.883, abs=1e-3)


def test_curve_returns_fractions(client):
    r = client.post("/api/curve", json={"totals": {"0006": 0.01}, "kinds": ["ABB"],
                                        "ph_min": 0, "ph_max": 14, "ph_points": 15})
    assert r.status_code == 200
    body = r.json()
    assert len(body["ph"]) == 15
    assert all(len(v) == 15 for v in body["fractions"].values())


def test_invalid_kind_rejected(client):
    r = client.post("/api/speciate", json={"totals": {"0002": 0.1}, "kinds": ["xyz"]})
    assert r.status_code == 422


def test_empty_totals_rejected(client):
    r = client.post("/api/speciate", json={"totals": {}})
    assert r.status_code == 422


def test_mc_endpoint(client):
    r = client.post("/api/mc", json={"totals": {"0002": 0.1, "0003": 0.1}, "kinds": ["ABB"],
                                     "sd": 0.02, "n": 30})
    assert r.status_code == 200
    assert r.json()["n_ok"] > 0
