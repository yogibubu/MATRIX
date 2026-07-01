from __future__ import annotations

import json
from pathlib import Path


def test_matrix_morpheus_imports_and_data_library():
    import matrix_morpheus as morpheus

    assert morpheus.SEMIEXP_JOB_SCHEMA == "oracle.semiexp.job.v1"
    assert morpheus.XYZIN_ISOTOPOLOGUES_SCHEMA == "oracle.xyz.isotopologues.v1"
    assert morpheus.DEFAULT_SE_GEOMETRY_LIBRARY.is_dir()
    assert (morpheus.DEFAULT_SE_GEOMETRY_LIBRARY / "manifest.csv").is_file()
    assert morpheus.DEFAULT_SEMIEXP_BENCHMARK_SNAPSHOT.parts[:2] == (
        "benchmarks",
        "semiexp_msr",
    )


def test_semiexp_paper_benchmark_paths_resolve_from_repo_root():
    from matrix_morpheus.paper_benchmarks import repository_root

    root = Path(__file__).resolve().parents[1]

    assert repository_root() == root


def test_semiexp_paper_benchmark_snapshot_generates_tables(tmp_path):
    from matrix_morpheus import generate_paper_benchmark_artifacts

    snapshot, artifacts = generate_paper_benchmark_artifacts(
        outdir=tmp_path / "paper",
        refresh_from_outputs=False,
    )

    assert snapshot["schema"] == "oracle.semiexp.paper_regression.v1"
    assert tuple(snapshot["cases"]) == (
        "glycolaldehyde",
        "glycine_I",
        "glycine_II",
        "cyclopentadiene",
        "nitrobenzene",
        "azulene",
        "norcamphor",
    )
    assert set(snapshot["planar_pair_diagnostics"]) == {"nitrobenzene", "azulene"}
    assert artifacts["summary_tex"].is_file()
    assert artifacts["planar_tex"].is_file()
    assert "Nitrobenzene" in artifacts["summary_tex"].read_text(encoding="utf-8")


def test_semiexp_job_accepts_oracle_and_legacy_schemas(tmp_path):
    from matrix_morpheus import read_semiexperimental_job

    body = """
title = "water inline"

[geometry]
units = "angstrom"
atoms = [
  ["O", 0.0, 0.0, 0.0],
  ["H", 0.0, 0.0, 0.9572],
  ["H", 0.9266, 0.0, -0.2396],
]
"""
    oracle_job = tmp_path / "water_oracle.mse.toml"
    legacy_job = tmp_path / "water_legacy.mse.toml"
    oracle_job.write_text('schema = "oracle.semiexp.job.v1"\n' + body, encoding="utf-8")
    legacy_job.write_text('schema = "merlino.semiexp.job.v1"\n' + body, encoding="utf-8")

    assert read_semiexperimental_job(oracle_job).geometry.source_format == "matrix_morpheus_job"
    assert read_semiexperimental_job(legacy_job).geometry.source_format == "matrix_morpheus_job"


def test_prepare_semiexp_xyzin_supports_custom_path(tmp_path):
    from matrix_morpheus import prepare_semiexperimental_xyzin, read_geometry_input

    root = Path(__file__).resolve().parents[1]
    xyzin = tmp_path / "water_container"
    result = prepare_semiexperimental_xyzin(
        root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz",
        observations_source=root
        / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml",
        xyzin_path=xyzin,
    )

    assert result.xyzin == xyzin
    assert "SCHEMA oracle.xyz.isotopologues.v1" in xyzin.read_text(encoding="utf-8")
    assert read_geometry_input(xyzin).source_format == "xyzin"


def test_oracle_semiexp_cli_runs_water_gic(tmp_path):
    from matrix_core.cli import main
    from matrix_morpheus import read_morpheus_section

    root = Path(__file__).resolve().parents[1]
    outdir = tmp_path / "run"
    xyzin = tmp_path / "xyzin"
    status = main(
        [
            "semiexp",
            "--xyz",
            str(root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz"),
            "--observations",
            str(root / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml"),
            "--xyzin",
            str(xyzin),
            "--outdir",
            str(outdir),
            "--coordinate-model",
            "gic",
            "--max-iter",
            "2",
        ],
        repo_root=root,
    )

    manifest = json.loads((outdir / "semiexp_manifest.json").read_text(encoding="utf-8"))
    checkpoint = json.loads((outdir / "semiexp_checkpoint.json").read_text(encoding="utf-8"))
    section = read_morpheus_section(xyzin)

    assert status == 0
    assert manifest["schema_version"] == "oracle.run.v1"
    assert checkpoint["schema"] == "oracle.semiexp.checkpoint.v1"
    assert section.status == "complete"
    assert section.run_dir == outdir
    assert section.manifest_path == outdir / "semiexp_manifest.json"
    assert section.html_report_path == outdir / "semiexp_report.html"
    assert section.latex_tables_path == outdir / "semiexp_tables.tex"
    assert section.coordinate_model == "gic"
    assert section.iterations <= 2
    assert section.parameter_count > 0
    assert "ORACLE semiexperimental equilibrium-geometry fit" in (
        outdir / "semiexp_report.txt"
    ).read_text(encoding="utf-8")
    assert (outdir / "semiexp_report.html").is_file()
