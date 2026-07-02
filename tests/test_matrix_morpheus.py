from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


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


def test_semiexp_covariance_uses_declared_weights_not_residual_rescaling():
    from matrix_morpheus.fit import _covariance

    jac = np.eye(2)
    small_residual = np.array([0.1, 0.2])
    large_residual = np.array([100.0, -50.0])

    expected = np.eye(2)
    assert np.allclose(_covariance(jac, small_residual), expected)
    assert np.allclose(_covariance(jac, large_residual), expected)


def test_semiexp_geometry_report_includes_neo_primitives_outside_covalent_graph():
    from matrix_morpheus.fit import _geometry_parameters
    from matrix_neo.survibfit.primitives import Primitive

    atoms = ("H", "O", "H")
    coords = np.array(
        [
            [-0.75, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.75, 0.0, 0.0],
        ],
        dtype=float,
    )

    rows = _geometry_parameters(atoms, coords, fit_prims=(Primitive("bond", (0, 2)),))

    assert any(row.kind == "neo_bond" and row.label == "NEO:R(1,3)" for row in rows)


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


def test_synthon_advisor_builds_topology_classes_for_equivalent_bonds():
    from matrix_morpheus import synthon_primitive_class_specs

    atoms = ("O", "H", "H")
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.9572],
            [0.9266, 0.0, -0.2396],
        ]
    )

    specs = synthon_primitive_class_specs(atoms, coords, include_angles=False)
    patterns = {pattern for spec in specs for pattern in spec.patterns}

    assert {"R(1,2)", "R(1,3)"} <= patterns


def test_synthon_advisor_refinement_level_changes_class_granularity():
    from matrix_morpheus import synthon_level_for_budget, synthon_primitive_class_specs

    atoms = ("C", "C", "O", "H", "H", "H", "H", "H")
    coords = np.array(
        [
            [0.000, 0.000, 0.000],
            [1.520, 0.000, 0.000],
            [2.900, 0.000, 0.000],
            [-0.540, 0.930, 0.000],
            [-0.540, -0.930, 0.000],
            [1.520, 1.030, 0.000],
            [1.520, -1.030, 0.000],
            [3.240, 0.900, 0.000],
        ]
    )

    coarse = synthon_primitive_class_specs(atoms, coords, level="coarse")
    fine = synthon_primitive_class_specs(atoms, coords, level="fine")

    assert len(coarse) <= len(fine)
    assert synthon_level_for_budget(4) == "coarse"
    assert synthon_level_for_budget(8) == "medium"
    assert synthon_level_for_budget(20) == "fine"


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


def test_morpheus_xyzin_config_reads_auto_predicates_and_synthon_advisor():
    from matrix_morpheus import parse_morpheus_input_config

    xyzin = """#MORPHEUS
COORDINATE_MODEL = gic
OBSERVABLE = moments
COMPONENTS = Ia,Ib,Ic
PREDICATES = INITIAL_GEOMETRY REFERENCE_LEVEL=high DISTANCE_SIGMA=AUTO ANGLE_SIGMA=AUTO
PREDICATE_SCOPE = xy_bonds,xh_bonds
PRIMITIVE_CLASS_ADVISOR = AUTO_SYNTHON INCLUDE=bonds MIN_GROUP_SIZE=2
SYNTHON_LEVEL = fine
PRIMITIVE_CLASS_BUDGET = auto
"""

    config = parse_morpheus_input_config(xyzin.splitlines()[1:])

    assert config.initial_geometry_predicates.distance_sigma_angstrom == pytest.approx(0.0015)
    assert config.initial_geometry_predicates.angle_sigma_degree == pytest.approx(0.15)
    assert config.synthon_primitive_classes.enabled
    assert config.synthon_primitive_classes.level == "fine"
    assert config.synthon_primitive_classes.include_bonds
    assert not config.synthon_primitive_classes.include_angles
    assert config.primitive_class_budget == "auto"


def test_morpheus_xyzin_config_normalizes_rotational_component_lists():
    from matrix_morpheus import parse_morpheus_input_config

    config = parse_morpheus_input_config(
        [
            "COORDINATE_MODEL = gic",
            "OBSERVABLE = rotational_constants",
            "COMPONENTS = A,B,C",
        ]
    )

    assert config.components == "ABC"


def test_initial_geometry_predicates_separate_xy_and_xh_bonds():
    from matrix_morpheus import initial_geometry_predicates

    atoms = ("C", "O", "H", "H")
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0],
            [-0.7, 0.9, 0.0],
            [1.8, 0.9, 0.0],
        ]
    )

    xy = initial_geometry_predicates(
        atoms,
        coords,
        distance_sigma_angstrom=0.003,
        angle_sigma_degree=0.3,
        dihedral_sigma_degree=0.5,
        scope=("xy_bonds",),
    )
    xh = initial_geometry_predicates(
        atoms,
        coords,
        distance_sigma_angstrom=0.003,
        angle_sigma_degree=0.3,
        dihedral_sigma_degree=0.5,
        scope=("xh_bonds",),
    )

    assert {item.label_pattern for item in xy} == {"R(1,2)"}
    assert {item.label_pattern for item in xh} == {"R(1,3)", "R(2,4)"}


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


def test_semiexp_sensitivity_advisor_generates_gic_predicates(tmp_path):
    from matrix_morpheus import (
        SemiexperimentalFitRequest,
        advise_semiexperimental_gic_sensitivity,
        read_observations_toml,
    )

    root = Path(__file__).resolve().parents[1]
    geometry = root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz"
    observations = read_observations_toml(
        root / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml"
    )
    request = SemiexperimentalFitRequest(geometry, observations)

    advisor = advise_semiexperimental_gic_sensitivity(
        request,
        fit_relative_threshold=1.1,
        fixed_relative_threshold=-1.0,
        min_fit_count=0,
    )
    csv_path = tmp_path / "advisor.csv"
    csv_path.write_text(advisor.csv, encoding="utf-8")

    assert advisor.predicate_count == len(advisor.predicates) > 0
    assert advisor.fit_count == 0
    assert all(item.source == "morpheus_sensitivity_advisor" for item in advisor.predicates)
    assert "relative_sensitivity,chemical_role,current_state,suggested_state" in (
        csv_path.read_text(encoding="utf-8")
    )


def test_semiexp_sensitivity_advisor_keeps_isotopologue_coverage():
    from matrix_morpheus import (
        SemiexperimentalFitRequest,
        advise_semiexperimental_gic_sensitivity,
        read_observations_toml,
    )

    root = Path(__file__).resolve().parents[1]
    geometry = root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz"
    observations = read_observations_toml(
        root / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml"
    )
    expanded = observations * 12
    request = SemiexperimentalFitRequest(geometry, expanded)

    advisor = advise_semiexperimental_gic_sensitivity(
        request,
        fit_relative_threshold=2.0,
    )

    assert advisor.fit_count == 2
    assert any(row.reason == "minimum_isotopologue_coverage" for row in advisor.rows)


def test_semiexp_sensitivity_advisor_prefers_soft_inter_coordinates():
    from matrix_morpheus.report import _sensitivity_selected_fit_labels

    candidates = [
        ("GIC001 NEO AStr0001 irrep=A R(1,2)", 1.0, 100.0, 1.0, ""),
        ("GIC002 NEO APsAn0001 irrep=A NONE", 0.0, 10.0, 0.1, ""),
        ("GIC003 NEO APsTo0001 irrep=A NONE", 0.0, 9.0, 0.09, ""),
        ("GIC004 NEO AStrD001 irrep=A LINEAR_COMBINATION", 0.0, 8.0, 0.08, ""),
        ("GIC005 NEO ABendD001 irrep=A LINEAR_COMBINATION", 0.0, 7.0, 0.07, ""),
    ]

    selected = _sensitivity_selected_fit_labels(
        candidates,
        free_budget=3,
        required_fit=3,
        fit_relative_threshold=0.5,
    )

    assert selected == {"GIC002", "GIC003", "GIC004"}


def test_semiexp_sensitivity_advisor_scales_soft_gic_predicates():
    from matrix_morpheus.report import (
        _sensitivity_fit_regularization_sigma,
        _sensitivity_predicate_sigma,
    )

    soft = "GIC001 NEO APsAn0001 irrep=A NONE"
    hard = "GIC002 NEO AStr0001 irrep=A R(1,2)"
    soft_sigma = _sensitivity_predicate_sigma(
        soft,
        relative=0.01,
        fit_relative_threshold=0.5,
        distance_sigma_angstrom=0.003,
        angle_sigma_degree=0.3,
        torsion_sigma_degree=0.5,
        soft_predicate_scale=0.5,
        null_predicate_scale=0.25,
    )
    hard_fit_sigma = _sensitivity_fit_regularization_sigma(
        hard,
        relative=1.0,
        fit_regularization_scale=2.0,
        distance_sigma_angstrom=0.003,
        angle_sigma_degree=0.3,
        torsion_sigma_degree=0.5,
    )
    soft_fit_sigma = _sensitivity_fit_regularization_sigma(
        soft,
        relative=1.0,
        fit_regularization_scale=2.0,
        distance_sigma_angstrom=0.003,
        angle_sigma_degree=0.3,
        torsion_sigma_degree=0.5,
    )

    assert soft_sigma == pytest.approx(np.deg2rad(0.5) * 0.5 * 0.25)
    assert hard_fit_sigma == 0.0
    assert soft_fit_sigma == pytest.approx(np.deg2rad(0.5) * 4.0)


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


def test_semiexp_cli_writes_sensitivity_advisor_csv(tmp_path):
    from matrix_core.cli import main

    root = Path(__file__).resolve().parents[1]
    outdir = tmp_path / "run"
    status = main(
        [
            "semiexp",
            "--xyz",
            str(root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz"),
            "--observations",
            str(root / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml"),
            "--outdir",
            str(outdir),
            "--coordinate-model",
            "gic",
            "--sensitivity-advisor",
            "--sensitivity-fit-threshold",
            "1.1",
            "--sensitivity-fixed-threshold",
            "-1.0",
            "--sensitivity-min-fit",
            "none",
            "--max-iter",
            "1",
        ],
        repo_root=root,
    )

    advisor_csv = (outdir / "semiexp_sensitivity_advisor.csv").read_text(encoding="utf-8")
    weight_diagnostics = (outdir / "semiexp_weight_diagnostics.csv").read_text(
        encoding="utf-8"
    )

    assert status == 0
    assert "below_fit_threshold" in advisor_csv
    assert ",predicate," in advisor_csv
    assert "morpheus_sensitivity_advisor" not in weight_diagnostics


def test_semiexp_cli_applies_sensitivity_advisor_only_when_requested(tmp_path):
    from matrix_core.cli import main

    root = Path(__file__).resolve().parents[1]
    outdir = tmp_path / "run"
    status = main(
        [
            "semiexp",
            "--xyz",
            str(root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz"),
            "--observations",
            str(root / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml"),
            "--outdir",
            str(outdir),
            "--coordinate-model",
            "gic",
            "--sensitivity-advisor",
            "--apply-sensitivity-advisor",
            "--sensitivity-fit-threshold",
            "1.1",
            "--sensitivity-fixed-threshold",
            "-1.0",
            "--sensitivity-min-fit",
            "none",
            "--max-iter",
            "1",
        ],
        repo_root=root,
    )

    weight_diagnostics = (outdir / "semiexp_weight_diagnostics.csv").read_text(
        encoding="utf-8"
    )

    assert status == 0
    assert "morpheus_sensitivity_advisor" in weight_diagnostics


def test_semiexp_sensitivity_gate_rejects_invalid_base_model(tmp_path):
    from matrix_core.cli import _sensitivity_safe_apply_gate

    def failing_fit(*_args, **_kwargs):
        raise RuntimeError("base model is underdetermined")

    gate = _sensitivity_safe_apply_gate(
        base_request=object(),
        candidate_request=object(),
        fit_semiexperimental_geometry=failing_fit,
        outdir=tmp_path,
        max_iter=1,
        step=1.0e-4,
        damping=1.0e-8,
        max_step=0.25,
        prune_condition=0.0,
        rot_rel_tol=0.02,
        rot_abs_tol=1.0e-3,
        condition_factor=10.0,
        max_bond_delta=0.01,
        max_angle_delta=1.0,
    )

    assert gate["accepted"] is False
    assert gate["reason"] == "base_model_invalid"


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
PREDICATE_SCOPE = xh_bonds
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
