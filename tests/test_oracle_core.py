from __future__ import annotations

from oracle_core import (
    ORACLE_MANIFEST_SCHEMA,
    build_run_manifest,
    ensure_workspace,
    replace_section,
    section_content,
)


def test_workspace_layout(tmp_path):
    layout = ensure_workspace(tmp_path / "project")
    for name in ("inputs", "runs", "outputs", "reports", "cache", "logs"):
        assert (layout.root / name).is_dir()


def test_manifest_schema(tmp_path):
    manifest = build_run_manifest(workflow="smoke", status="completed", run_dir=tmp_path)
    assert manifest.schema_version == ORACLE_MANIFEST_SCHEMA
    assert manifest.to_dict()["workflow"] == "smoke"


def test_section_replacement_preserves_other_sections(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0 0 0",
                "H 0 0 1",
                "H 0 1 0",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
                "BOND 1 2",
                "",
                "#GIC",
                "SCHEMA oracle.xyz.gic.v1",
                "OLD true",
                "",
                "#MORPHEUS",
                "SCHEMA oracle.xyz.morpheus.v1",
                "STATUS draft",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    replace_section(path, "GIC", ["SCHEMA oracle.xyz.gic.v1", "UPDATED true"])
    lines = path.read_text(encoding="utf-8").splitlines()

    assert section_content(lines, "TOPOLOGY") == [
        "SCHEMA oracle.xyz.topology.v1",
        "BOND 1 2",
        "",
    ]
    assert section_content(lines, "GIC") == [
        "SCHEMA oracle.xyz.gic.v1",
        "UPDATED true",
    ]
    assert section_content(lines, "MORPHEUS") == [
        "SCHEMA oracle.xyz.morpheus.v1",
        "STATUS draft",
    ]
