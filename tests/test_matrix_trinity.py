from __future__ import annotations

import json
from pathlib import Path


def test_trinity_section_roundtrip(tmp_path):
    from matrix_trinity import (
        ORACLE_XYZ_TRINITY_SCHEMA,
        TrinitySection,
        read_trinity_section,
        write_trinity_section,
    )

    xyzin = tmp_path / "molecule.xyzin"
    run_dir = tmp_path / "trinity"
    xyzin.write_text("1\nh\nH 0.0 0.0 0.0\n", encoding="utf-8")

    write_trinity_section(
        xyzin,
        TrinitySection(
            source_path=xyzin,
            run_dir=run_dir,
            manifest_path=run_dir / "trinity_manifest.json",
            engine_command="engine --gradient",
            coordinate_model="gic",
            active_space="total_symmetric",
            max_steps=12,
            trust_radius=0.15,
            gradient_tolerance=2.0e-5,
            step_tolerance=3.0e-5,
            energy_tolerance=4.0e-8,
            trajectory_path=run_dir / "traj.xyz",
            final_geometry_path=run_dir / "final.xyz",
            energy_gradient_log_path=run_dir / "eg.jsonl",
        ),
    )

    section = read_trinity_section(xyzin)

    assert section.schema == ORACLE_XYZ_TRINITY_SCHEMA
    assert section.source_path == xyzin
    assert section.run_dir == run_dir
    assert section.engine_command == "engine --gradient"
    assert section.coordinate_model == "gic"
    assert section.active_space == "total_symmetric"
    assert section.max_steps == 12
    assert section.trust_radius == 0.15
    assert section.gradient_tolerance == 2.0e-5
    assert section.final_geometry_path == run_dir / "final.xyz"


def test_trinity_cli_prepare_and_status(tmp_path, capsys):
    from matrix_core.cli import main
    from matrix_trinity import read_trinity_section

    root = Path(__file__).resolve().parents[1]
    xyzin = tmp_path / "molecule.xyzin"
    run_dir = tmp_path / "trinity_run"
    xyzin.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0.0 0.0 0.0",
                "",
                "#BASIC",
                "SCHEMA oracle.xyz.basic.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rc = main(
        [
            "trinity",
            "prepare",
            str(xyzin),
            "--run-dir",
            str(run_dir),
            "--engine-command",
            "external-engine --input step.xyz --gradient gradient.json",
            "--max-steps",
            "7",
            "--trust-radius",
            "0.12",
        ],
        repo_root=root,
    )
    out = capsys.readouterr().out
    section = read_trinity_section(xyzin)
    manifest = json.loads((run_dir / "trinity_manifest.json").read_text(encoding="utf-8"))

    assert rc == 0
    assert "Updated #TRINITY" in out
    assert section.status == "prepared"
    assert section.run_dir == run_dir
    assert section.manifest_path == run_dir / "trinity_manifest.json"
    assert section.engine_command.startswith("external-engine")
    assert section.max_steps == 7
    assert section.trust_radius == 0.12
    assert manifest["workflow"] == "trinity"
    assert manifest["status"] == "prepared"
    assert manifest["parameters"]["max_steps"] == 7
    assert manifest["backend"]["engine_command"].startswith("external-engine")

    rc = main(["trinity", "status", str(xyzin)], repo_root=root)
    status_out = capsys.readouterr().out

    assert rc == 0
    assert "coordinate model: gic" in status_out
    assert "max steps: 7" in status_out
