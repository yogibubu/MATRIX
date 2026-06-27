from pathlib import Path

import numpy as np

from matrix_chem import preprocess_to_enriched_xyz, read_enriched_xyz, write_validation_section
from matrix_gaussian import (
    check_gaussian_readallgic_log,
    hessian_input_from_gaussian_log,
    promote_gaussian_log_hessian_to_xyzin,
    read_gaussian_log_geometry,
    read_gaussian_log_cartesian_hessian,
    read_gaussian_log_normal_modes,
    summarize_gaussian_log,
)
from matrix_gf import run_xyzin_gf_report_from_xyzin
from matrix_gf.harmonic import HESSIAN_EIGENVALUE_TO_CM
from matrix_neo import write_gicforge_build_sections
from matrix_qm import read_normal_modes_section


ROOT = Path(__file__).resolve().parents[1]
PYRROLE = ROOT / "examples" / "gaussian_readallgic" / "pyrrole"


def test_pyrrole_readallgic_example_log_is_valid() -> None:
    check = check_gaussian_readallgic_log(
        PYRROLE / "pyrrole.log",
        expected_point_group="C2V",
        expected_rank=24,
        expected_frozen_labels=(
            "B2Str001",
            "B2CyBe001",
            "A2RPck001",
            "B1RPck001",
            "A2OuPl001",
            "B1OuPl003",
        ),
        expected_active_labels=("A1Str001", "A1Bend001", "A1CyBe001"),
        require_frequency=True,
        require_no_imaginary=True,
    )

    assert check.ok, check.errors
    assert check.normal_termination_count >= 2
    assert check.route_has_readallgic
    assert check.optimization_completed
    assert check.stationary_point
    assert check.frequency_count == 24
    assert check.n_imag == 0


def test_pyrrole_readallgic_input_freezes_non_total_coordinates() -> None:
    text = (PYRROLE / "pyrrole.gjf").read_text(encoding="utf-8").lower()

    assert "opt=readallgic" in text
    assert "geom=readallgic" not in text
    assert "gicallsym" not in text
    assert "a1str001(frozen)" not in text
    assert "a1bend001(frozen)" not in text
    assert "a1cybe001(frozen)" not in text
    assert "a2rpck001(frozen)" in text
    assert "b1rpck001(frozen)" in text
    assert "a2oupl001(frozen)" in text
    assert "b1oupl003(frozen)" in text


def test_pyrrole_log_hessian_reconstructs_gaussian_frequencies() -> None:
    raw_hessian = read_gaussian_log_cartesian_hessian(PYRROLE / "pyrrole.log")
    hessian_input = hessian_input_from_gaussian_log(PYRROLE / "pyrrole.log")
    weights = 1.0 / np.sqrt(np.repeat(hessian_input.masses_amu, 3))
    mass_weighted = hessian_input.cartesian_hessian * weights[:, None] * weights[None, :]
    eigenvalues = np.linalg.eigvalsh(0.5 * (mass_weighted + mass_weighted.T))
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * HESSIAN_EIGENVALUE_TO_CM
    gaussian_frequencies = np.asarray(
        summarize_gaussian_log(PYRROLE / "pyrrole.log").frequencies_cm
    )

    assert hessian_input.source == "gaussian-log"
    assert raw_hessian.shape == hessian_input.cartesian_hessian.shape
    assert np.allclose(raw_hessian, hessian_input.cartesian_hessian)
    assert hessian_input.cartesian_hessian.shape == (30, 30)
    assert np.allclose(frequencies[6:], gaussian_frequencies, atol=0.4)


def test_pyrrole_log_geometry_is_final_archive_geometry() -> None:
    geometry = read_gaussian_log_geometry(PYRROLE / "pyrrole.log")
    hessian_input = hessian_input_from_gaussian_log(PYRROLE / "pyrrole.log")

    assert geometry.source_format == "gaussian_log_archive"
    assert geometry.metadata["orientation"] == "archive_original_axes"
    assert np.allclose(
        geometry.coordinates_angstrom,
        hessian_input.cartesian_coordinates_bohr * 0.52917721092,
    )


def test_pyrrole_gf_uses_log_hessian_and_symmetrized_gics(tmp_path: Path) -> None:
    xyzin = tmp_path / "pyrrole.xyzin"
    preprocess_to_enriched_xyz(PYRROLE / "pyrrole.log", xyzin, source_kind="gaussian")
    write_validation_section(xyzin)
    write_gicforge_build_sections(xyzin, symmetrize=True)
    promotion = promote_gaussian_log_hessian_to_xyzin(PYRROLE / "pyrrole.log", xyzin)

    report = run_xyzin_gf_report_from_xyzin(xyzin, block_by_irrep=True)
    gaussian_frequencies = np.asarray(
        summarize_gaussian_log(PYRROLE / "pyrrole.log").frequencies_cm
    )
    normal_modes = read_normal_modes_section(xyzin)
    geometry = read_enriched_xyz(xyzin)
    ped = report.result.ped.values

    assert promotion.wrote_cartesian_hessian is True
    assert promotion.wrote_normal_modes is True
    assert geometry.source_format == "enriched_xyz"
    assert report.result.point_group == "C2v"
    assert report.result.symmetrized_gics is True
    assert report.result.matrix_model == "FULL+IRREP_BLOCKS"
    assert report.result.gic_names[:3] == ("A1Str001", "A1Str002", "A1Str003")
    assert "A2RPck001" in report.result.gic_names
    assert "B1OuPl003" in report.result.gic_names
    assert len(report.result.block_labels) == len(report.result.frequencies_cm)
    assert normal_modes.modes.shape == (24, 30)
    assert np.allclose(normal_modes.frequencies_cm, gaussian_frequencies)
    assert report.frequency_comparison is not None
    assert report.geometry_comparison is not None
    assert "Frequency check vs" in report.text
    assert "Geometry check:" in report.text
    assert report.geometry_comparison.raw_max_angstrom < 1.0e-7
    assert report.geometry_comparison.aligned_max_angstrom < 1.0e-7
    assert report.frequency_comparison.max_abs_delta_cm < 0.02
    assert np.allclose(
        np.sort(report.result.frequencies_cm),
        np.sort(gaussian_frequencies),
        atol=0.02,
    )
    assert np.all(ped >= -1.0e-12)
    assert np.allclose(np.sum(ped, axis=0), 100.0, atol=1.0e-8)
    gic_irreps = np.asarray(report.result.gic_irreps)
    for mode, irrep in enumerate(report.result.block_labels):
        off_block_ped = ped[gic_irreps != irrep, mode]
        assert np.max(off_block_ped, initial=0.0) < 1.0e-8


def test_pyrrole_log_normal_modes_parser_reads_gaussian_table() -> None:
    normal_modes = read_gaussian_log_normal_modes(PYRROLE / "pyrrole.log")
    gaussian_frequencies = np.asarray(
        summarize_gaussian_log(PYRROLE / "pyrrole.log").frequencies_cm
    )

    assert normal_modes.modes.shape == (24, 30)
    assert normal_modes.rotated_to_archive_axes is True
    assert np.allclose(normal_modes.frequencies_cm, gaussian_frequencies)
