from __future__ import annotations

import json
from pathlib import Path

from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
from matrix_fragments import write_interaction_center_section
from matrix_link import rdkit_available
from matrix_neo import (
    DEFAULT_FORTRAN_AUDIT_MOLECULES,
    audit_gic_corpus_geometry,
    discover_gic_corpus,
    summarize_gic_corpus,
    write_gicforge_build_sections,
)


CORPUS = Path(__file__).resolve().parent / "fixtures" / "test_molecules" / "molecules"
GOLDEN_CORPUS = (
    Path(__file__).resolve().parent / "fixtures" / "golden_corpus" / "neo_gic_golden.json"
)
GOLDEN_SALC_COEFFICIENTS = (
    Path(__file__).resolve().parent
    / "fixtures"
    / "golden_corpus"
    / "neo_gic_salc_coefficients.json"
)
ENV_HELPERS = Path(__file__).resolve().parents[1] / "scripts" / "matrix_env.sh"


def test_gic_regression_corpus_is_available():
    inputs = sorted(CORPUS.glob("*.inp"))
    names = {path.name for path in inputs}

    assert len(inputs) >= 100
    assert {
        "benzene.inp",
        "cyclopentadiene_B3.inp",
        "nitrobenzene.inp",
        "azulene.inp",
        "norbornane.inp",
    } <= names


def test_official_neo_gic_golden_corpus_registry_is_stable():
    root = Path(__file__).resolve().parents[1]
    registry = json.loads(GOLDEN_CORPUS.read_text(encoding="utf-8"))
    entries = registry["entries"]
    ids = {entry["id"] for entry in entries}
    roles = {role for entry in entries for role in entry["roles"]}

    assert registry["schema"] == "matrix.golden_corpus.v1"
    assert len(entries) >= 20
    assert len(ids) == len(entries)
    assert {
        "benzene",
        "pyrrole",
        "azulene",
        "norbornane",
        "norbornene",
        "norbornadiene",
        "norcamphor",
        "thujone",
        "ribose",
        "cubane",
        "ferrocene_d5h",
        "ferrocene_d5d",
        "formic_acid_water_hbond",
        "spiro",
        "pyrene",
    } <= ids
    assert {
        "ring",
        "fused_ring",
        "polycyclic",
        "bridged_ring",
        "spiro_ring",
        "hbond",
        "pseudo_bond",
        "metal_center",
        "special_coordinates",
        "python_fortran_parity",
    } <= roles
    for entry in entries:
        assert (root / entry["path"]).is_file(), entry["path"]


def test_fortran_audit_default_covers_golden_gic_closure_roles():
    registry = json.loads(GOLDEN_CORPUS.read_text(encoding="utf-8"))
    closure_roles = {
        "fused_ring",
        "bridged_ring",
        "spiro_ring",
        "python_fortran_parity",
        "metal_center",
        "ring_center",
        "d5h_symmetry",
        "d5d_symmetry",
    }
    required = {
        Path(entry["path"]).name
        for entry in registry["entries"]
        if closure_roles.intersection(entry["roles"])
        and Path(entry["path"]).parent.name == "molecules"
    }

    assert required <= set(DEFAULT_FORTRAN_AUDIT_MOLECULES)


def test_neo_gic_salc_coefficient_snapshots_are_stable(tmp_path):
    root = Path(__file__).resolve().parents[1]
    snapshot = json.loads(GOLDEN_SALC_COEFFICIENTS.read_text(encoding="utf-8"))

    assert snapshot["schema"] == "matrix.neo.gic_salc_coefficients.v1"
    assert snapshot["rounding_decimals"] == 12
    assert len(snapshot["entries"]) >= 8

    for entry in snapshot["entries"]:
        xyzin = tmp_path / f"{entry['id']}.xyzin"
        preprocess_to_enriched_xyz(root / entry["source"], xyzin)
        write_validation_section(xyzin)
        if str(entry["id"]).startswith("ferrocene"):
            write_interaction_center_section(xyzin)
        definition = write_gicforge_build_sections(xyzin, symmetrize=True)

        current = [
            {
                "name": gic.name,
                "family": gic.family,
                "irrep": gic.irrep,
                "coefficients": [
                    [primitive_id, round(float(coefficient), snapshot["rounding_decimals"])]
                    for primitive_id, coefficient in gic.coefficients
                ],
            }
            for gic in definition.gics
            if len(gic.coefficients) > 1
        ]

        assert definition.point_group == entry["point_group"], entry["id"]
        assert definition.rank == entry["rank"], entry["id"]
        assert definition.target_rank == entry["target_rank"], entry["id"]
        assert definition.symmetry_diagnostics is not None
        assert definition.symmetry_diagnostics.method == entry["symmetry_method"], entry["id"]
        assert len(current) == entry["salc_count"], entry["id"]
        assert current == entry["salcs"], entry["id"]


def test_gic_regression_corpus_keeps_qm_adapter_outputs():
    assert (CORPUS / "apinene.log").is_file()
    assert (CORPUS / "apinene.out").is_file()
    assert (CORPUS / "c6h5.fchk").is_file()
    assert (CORPUS / "c6h5.gjf").is_file()


def test_gic_regression_corpus_inventory_classifies_files():
    summary = summarize_gic_corpus(CORPUS)
    inp_entries = discover_gic_corpus(CORPUS, suffixes=["inp"])

    assert summary.total_files == 156
    assert summary.suffix_counts[".inp"] == 129
    assert summary.role_counts["legacy_gic_input"] == 129
    assert len(inp_entries) == 129
    assert {entry.role for entry in inp_entries} == {"legacy_gic_input"}


def test_gic_regression_corpus_geometry_audit_tracks_parser_budget():
    audit = audit_gic_corpus_geometry(CORPUS)
    failures = {entry.name for entry in audit.entries if not entry.passed}

    assert audit.total_files == 132
    assert audit.passed_files >= 115
    assert audit.failed_files <= 15
    assert audit.source_format_counts["gaussian_cartesian_input"] == 102
    if rdkit_available():
        assert audit.source_format_counts["gaussian_zmatrix_input"] <= 15
    else:
        assert audit.source_format_counts["gaussian_zmatrix_input"] == 15
    if rdkit_available():
        assert "testvib.inp" in failures
    else:
        assert audit.passed_files == 117
        assert audit.failed_files == 15
        assert audit.error_counts == {"GeometryParseError": 1, "RDKitUnavailableError": 14}
        assert {"azulene.inp", "pyrrole_smile1.inp", "testvib.inp"} <= failures


def test_matrix_environment_helpers_define_primary_and_legacy_commands():
    text = ENV_HELPERS.read_text(encoding="utf-8")

    for name in (
        "matrix-set()",
        "matrix-unset()",
        "matrix-run()",
        "matrix-run-bg()",
        "matrix-run-check()",
        "matrix-test-all()",
        "matrix-clean()",
        "matrix-create-venv()",
        "matrix-save-shell-state()",
        "matrix-install-runtime-deps()",
        "matrix-gic-corpus-list()",
        "matrix-gic-corpus-summary()",
        "matrix-gic-corpus-audit()",
        "oracle-set()",
        "oracle-unset()",
        "oracle-run()",
        "oracle-run-bg()",
        "oracle-run-check()",
        "oracle-test-all()",
        "oracle-clean()",
        "oracle-create-venv()",
        "oracle-save-shell-state()",
        "oracle-install-runtime-deps()",
        "oracle-gic-corpus-list()",
        "oracle-gic-corpus-summary()",
        "oracle-gic-corpus-audit()",
    ):
        assert name in text
    assert "MATRIX_AUTO_CREATE_VENV" in text
    assert "MATRIX_AUTO_INSTALL_RUNTIME_DEPS" in text
    assert "MATRIX_ENV_ACTIVE" in text
    assert "rdkit" in text
