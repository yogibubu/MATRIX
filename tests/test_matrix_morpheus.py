from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np


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
    from matrix_morpheus import generate_paper_benchmark_artifacts, validate_paper_run_outputs

    root = Path(__file__).resolve().parents[1]
    snapshot, artifacts = generate_paper_benchmark_artifacts(
        outdir=tmp_path / "paper",
        refresh_from_outputs=False,
    )

    assert snapshot["schema"] == "oracle.semiexp.paper_regression.v1"
    assert tuple(snapshot["cases"]) == (
        "glycolaldehyde",
        "cyclopentadiene",
        "nitrobenzene",
        "p-EBN",
        "azulene",
        "norcamphor",
        "glycine_I",
        "glycine_II",
    )
    assert set(snapshot["planar_pair_diagnostics"]) == {"nitrobenzene", "azulene"}
    assert artifacts["summary_tex"].is_file()
    assert artifacts["planar_tex"].is_file()
    run_outputs = validate_paper_run_outputs(snapshot)
    assert set(run_outputs) == {"complete", "missing", "incomplete"}
    assert not run_outputs["incomplete"]
    assert "Nitrobenzene" in artifacts["summary_tex"].read_text(encoding="utf-8")
    paper_generated = root / "docs/papers/morpheus_semiexp/generated"
    for name, artifact in (
        ("paper_benchmark_summary.csv", artifacts["summary_csv"]),
        ("paper_planar_pair_diagnostics.csv", artifacts["planar_csv"]),
        ("benchmark_summary.tex", artifacts["summary_tex"]),
        ("planar_pair_diagnostics.tex", artifacts["planar_tex"]),
    ):
        assert artifact.read_text(encoding="utf-8") == (
            paper_generated / name
        ).read_text(encoding="utf-8")


def test_semiexp_weight_diagnostics_track_robust_and_influence_terms():
    from matrix_morpheus.statistics import weight_diagnostic_rows, weight_diagnostics_csv

    model = SimpleNamespace(
        labels=(("parent", "A"), ("d1", "A"), ("predicate", "CC_class")),
        observed=np.array([10.0, 20.0, 0.0]),
        weights=np.array([1.0, 4.0, 9.0]),
        n_experimental_rows=2,
    )
    calculated = np.array([9.0, 18.0, 0.5])
    residual = model.observed - calculated
    robust_sqrt = np.array([1.0, 0.5, 1.0])
    weighted_residual = residual * np.sqrt(model.weights) * robust_sqrt
    rows = weight_diagnostic_rows(
        model,
        calculated,
        residual,
        np.eye(3, 2),
        weighted_residual,
        robust_sqrt,
    )

    assert len(rows) == 3
    assert rows[1].sigma == 0.5
    assert rows[1].robust_weight == 0.25
    assert rows[1].effective_weight == 1.0
    assert rows[2].kind == "predicate"
    assert rows[0].leverage > 0.0
    assert "total_weight_fraction" in weight_diagnostics_csv(rows).splitlines()[0]


def test_semiexp_rejects_underdetermined_free_refinement():
    import pytest

    from matrix_core import ScientificValidationError
    from matrix_morpheus.fit import MeasurementModel, _validate_observation_budget

    model = MeasurementModel(
        observable="moments",
        components=("Ia", "Ib", "Ic"),
        labels=(("parent", "Ia"), ("parent", "Ib"), ("parent", "Ic")),
        observed=np.ones(3),
        weights=np.ones(3),
        n_experimental_rows=3,
        planar=False,
    )

    with pytest.raises(ScientificValidationError, match="Underdetermined MORPHEUS"):
        _validate_observation_budget(model, 4, coordinate_model="gic")


def test_legacy_msr_zmatrix_constraints_import_for_p_ebn():
    from matrix_morpheus.msr_legacy import read_msr_legacy_input

    root = Path(__file__).resolve().parents[1]
    legacy = read_msr_legacy_input(
        root / "benchmarks/semiexp_msr/inputs/legacy_msr_import/p-EBN.msr.inp"
    )

    assert legacy.geometry.source_format == "msr_legacy_zmatrix"
    assert len(legacy.geometry.atoms) == 15
    assert len(legacy.observations) == 11
    assert "CC2=R(1,2)" in legacy.geometry.fixed_parameters
    assert "cc2 - cc4 Frozen" in legacy.geometry.fixed_parameters


def test_legacy_msr_parser_accepts_dvib_alias(tmp_path):
    from matrix_morpheus.msr_legacy import read_msr_legacy_input

    root = Path(__file__).resolve().parents[1]
    source = root / "benchmarks/semiexp_msr/inputs/legacy_msr_import/p-EBN.msr.inp"
    target = tmp_path / "p-EBN_dvib.msr.inp"
    target.write_text(source.read_text(encoding="utf-8").replace("dbvib", "dvib"), encoding="utf-8")

    legacy = read_msr_legacy_input(target)

    assert len(legacy.observations) == 11
    assert legacy.observations[0].correction.source == "MSR dbvib"
    assert any(value != 0.0 for value in legacy.observations[0].correction.as_tuple())


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


def test_primitive_class_advisor_maps_user_classes_to_gics():
    from matrix_morpheus import derive_primitive_class_plan, parse_primitive_class_spec

    labels = (
        "GIC001 AStr [ 0.75*R(1,2)+0.10*R(7,13)]",
        "GIC002 AStr [ 0.25*R(1,2)+0.74*R(7,13)]",
        "GIC003 AStr [ 0.74*R(1,2)+0.74*R(7,13)]",
        "GIC004 AStr [ 0.60*R(1,2)]",
    )
    classes = (
        parse_primitive_class_spec("CC_skeleton:R(2,1)"),
        parse_primitive_class_spec("CO_stretches:R(7,13)"),
    )

    plan = derive_primitive_class_plan(
        labels,
        classes,
        min_fraction=0.70,
        cross_fraction_max=0.20,
    )

    assert [(item.name, item.patterns) for item in plan.parameter_classes] == [
        ("CO_stretches", ("GIC002", "GIC003")),
        ("CC_skeleton", ("GIC001",)),
    ]
    assert ("CC_skeleton", 1, 0.75) in plan.class_support
    assert "GIC003" in plan.ambiguous_labels
    assert plan.fixed_patterns == ("GIC004",)
    assert plan.rejected_labels == ("GIC004",)


def test_primitive_class_advisor_respects_data_budget():
    from matrix_morpheus import derive_primitive_class_plan, parse_primitive_class_spec

    labels = (
        "GIC001 AStr [ 0.90*R(1,2)]",
        "GIC002 AStr [ 0.80*R(2,3)]",
        "GIC003 AStr [ 0.75*R(3,4)]",
        "GIC004 AStr [ 0.78*R(4,5)]",
    )
    classes = (
        parse_primitive_class_spec("class_a:R(1,2)|R(2,3)"),
        parse_primitive_class_spec("class_b:R(3,4)"),
        parse_primitive_class_spec("class_c:R(4,5)"),
    )

    plan = derive_primitive_class_plan(labels, classes, max_classes=1)

    assert [(item.name, item.patterns) for item in plan.parameter_classes] == [
        ("class_a", ("GIC001", "GIC002")),
    ]
    assert plan.budget_limited_classes == ("class_b", "class_c")
    assert plan.fixed_patterns == ("GIC003", "GIC004")


def test_primitive_class_advisor_supports_nonbond_primitives():
    from matrix_morpheus import derive_primitive_class_plan, parse_primitive_class_spec

    labels = (
        "GIC001 ABend [ 0.80*A(1,2,3)]",
        "GIC002 ATor [ 0.72*D(1,2,3,4)]",
        "GIC003 ABend [ 0.69*A(5,6,7)]",
    )
    classes = (
        parse_primitive_class_spec("bends:A(3,2,1)"),
        parse_primitive_class_spec("torsions:D(4,3,2,1)"),
    )

    plan = derive_primitive_class_plan(labels, classes)

    assert [(item.name, item.patterns) for item in plan.parameter_classes] == [
        ("bends", ("GIC001",)),
        ("torsions", ("GIC002",)),
    ]
    assert plan.fixed_patterns == ("GIC003",)


def test_semiexp_advisor_suggests_disjoint_dominant_gic_classes():
    from matrix_morpheus import suggest_parameter_classes

    labels = (
        "GIC001 AStr [ 0.90*R(1,2)+0.10*R(2,3)]",
        "GIC002 AStr [ 0.88*R(2,3)+0.05*R(1,4)]",
        "GIC003 AAng [ 0.82*A(1,2,3)]",
        "GIC004 AAng [ 0.78*A(2,3,4)]",
        "GIC005 AAng [ 0.91*A(1,5,2)]",
    )

    classes = suggest_parameter_classes(("C", "C", "C", "C", "H"), labels)

    assert (
        "CC_stretches",
        ("GIC001", "GIC002"),
        "shared",
    ) in [(item.name, item.patterns, item.mode) for item in classes]
    assert (
        "CCC_bends",
        ("GIC003", "GIC004"),
        "shared",
    ) in [(item.name, item.patterns, item.mode) for item in classes]
    assert (
        "XH_angle_directions",
        ("GIC005",),
        "fixed",
    ) in [(item.name, item.patterns, item.mode) for item in classes]


def test_kraitchman_seed_predicates_can_use_partial_seed_atoms():
    from matrix_morpheus import (
        kraitchman_seed_predicates,
        read_geometry_input,
        read_observations_toml,
    )

    root = Path(__file__).resolve().parents[1]
    geometry = read_geometry_input(
        root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz"
    )
    observations = read_observations_toml(
        root / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml"
    )

    conservative = kraitchman_seed_predicates(
        geometry.atoms,
        geometry.coordinates_angstrom,
        observations,
    )
    partial = kraitchman_seed_predicates(
        geometry.atoms,
        geometry.coordinates_angstrom,
        observations,
        require_all_atoms_seeded=False,
    )

    assert conservative == ()
    assert {item.label_pattern for item in partial} == {"R(1,2)", "R(1,3)", "A(2,1,3)"}
    assert all(item.source.startswith("kraitchman_seed") for item in partial)


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
    weight_diagnostics = (outdir / "semiexp_weight_diagnostics.csv").read_text(
        encoding="utf-8"
    )

    assert status == 0
    assert manifest["schema_version"] == "oracle.run.v1"
    assert manifest["outputs"]["weight_diagnostics"].endswith(
        "semiexp_weight_diagnostics.csv"
    )
    assert checkpoint["schema"] == "oracle.semiexp.checkpoint.v1"
    assert section.status == "complete"
    assert section.run_dir == outdir
    assert section.manifest_path == outdir / "semiexp_manifest.json"
    assert section.html_report_path == outdir / "semiexp_report.html"
    assert section.latex_tables_path == outdir / "semiexp_tables.tex"
    assert section.coordinate_model == "gic"
    assert section.iterations <= 2
    assert section.parameter_count > 0
    assert "MATRIX/MORPHEUS semiexperimental equilibrium-geometry fit" in (
        outdir / "semiexp_report.txt"
    ).read_text(encoding="utf-8")
    assert "Publishable Summary" in (outdir / "semiexp_report.html").read_text(
        encoding="utf-8"
    )
    assert (outdir / "semiexp_report.html").is_file()
    assert "base_weight,robust_weight,effective_weight" in weight_diagnostics
    assert "studentized_residual,cooks_distance" in weight_diagnostics


def test_semiexp_cli_accepts_standard_xyzin_morpheus_input(tmp_path):
    from matrix_core.cli import main
    from matrix_morpheus import prepare_semiexperimental_xyzin, read_morpheus_section

    root = Path(__file__).resolve().parents[1]
    xyzin = tmp_path / "water.xyzin"
    prepare_semiexperimental_xyzin(
        root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz",
        observations_source=root
        / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml",
        xyzin_path=xyzin,
    )
    xyzin.write_text(
        xyzin.read_text(encoding="utf-8")
        + """

#MORPHEUS
COORDINATE_MODEL = cartesian_symmetry
OBSERVABLE = moments
COMPONENTS = Ia,Ib
FIT_COORDINATES = cartesian_symmetry
PREDICATES = INITIAL_GEOMETRY DISTANCE_SIGMA=0.003 ANGLE_SIGMA=0.3 DIHEDRAL_SIGMA=0.5
PREDICATE_SCOPE = oh_bonds
""",
        encoding="utf-8",
    )
    outdir = tmp_path / "run"

    status = main(
        [
            "semiexp",
            "--xyzin",
            str(xyzin),
            "--outdir",
            str(outdir),
            "--max-iter",
            "1",
            "--no-write-section",
        ],
        repo_root=root,
    )

    report = (outdir / "semiexp_report.txt").read_text(encoding="utf-8")
    weights = (outdir / "semiexp_weight_diagnostics.csv").read_text(encoding="utf-8")
    section = read_morpheus_section(xyzin)

    assert status == 0
    assert "coordinate_model = cartesian_symmetry" in report
    assert "qm_predicate = R(1,2); value=0.9572; sigma=0.003" in report
    assert "qm_predicate = R(1,3); value=" in report
    assert ",predicate," in weights
    assert section.coordinate_model == "cartesian_symmetry"
