from __future__ import annotations

from oracle_core import ORACLE_MANIFEST_SCHEMA, build_run_manifest, ensure_workspace


def test_workspace_layout(tmp_path):
    layout = ensure_workspace(tmp_path / "project")
    for name in ("inputs", "runs", "outputs", "reports", "cache", "logs"):
        assert (layout.root / name).is_dir()


def test_manifest_schema(tmp_path):
    manifest = build_run_manifest(workflow="smoke", status="completed", run_dir=tmp_path)
    assert manifest.schema_version == ORACLE_MANIFEST_SCHEMA
    assert manifest.to_dict()["workflow"] == "smoke"

