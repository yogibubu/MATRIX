from __future__ import annotations

from pathlib import Path

import numpy as np

from oracle_chem import preprocess_to_enriched_xyz, write_validation_section
from oracle_gaussian import hessian_input_from_gaussian_fchk, lower_to_symmetric, read_gaussian_fchk
from oracle_gf import run_xyzin_gf_report_from_fchk, solve_wilson_gf, write_csv_tables
from oracle_gicforge import write_gicforge_build_sections


ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ROOT / "tests" / "fixtures" / "test_molecules" / "molecules"
GF_FIXTURES = ROOT / "tests" / "fixtures" / "gf"


def test_gaussian_fchk_hessian_adapter_reads_merlino_blocks():
    data = read_gaussian_fchk(GF_FIXTURES / "h2o.fchk")
    hessian = lower_to_symmetric(data.cartesian_hessian_lower)
    canonical = hessian_input_from_gaussian_fchk(GF_FIXTURES / "h2o.fchk")

    assert data.atomic_numbers.tolist() == [1, 8, 1]
    assert data.cartesian_coordinates_bohr.shape == (3, 3)
    assert hessian.shape == (9, 9)
    assert np.allclose(hessian, hessian.T)
    assert canonical.cartesian_hessian.shape == (9, 9)


def test_wilson_gf_solver_matches_diagonal_reference():
    result = solve_wilson_gf(
        np.diag([4.0, 9.0]),
        np.eye(2),
        scale_to_cm=False,
    )

    assert result.eigenvalues.tolist() == [4.0, 9.0]
    assert result.frequencies_cm.tolist() == [2.0, 3.0]


def test_xyzin_gf_report_runs_from_fchk_and_frozen_gics(tmp_path):
    source = MOLECULES / "h2ocart.inp"
    xyzin = tmp_path / "h2o.xyzin"
    scale = tmp_path / "scale.txt"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    scale.write_text("GIC003 0.90\n", encoding="utf-8")

    report = run_xyzin_gf_report_from_fchk(
        GF_FIXTURES / "h2o.fchk",
        xyzin,
        scale_path=scale,
    )
    written = write_csv_tables(report, tmp_path / "csv", prefix="gic_gf")

    assert "GF/PED from ORACLE non-redundant GICs" in report.text
    assert "Frozen xyzin:" in report.text
    assert "Pulay Hessian scaling: applied" in report.text
    assert len(report.result.frequencies_cm) == definition.rank
    assert np.all(np.isfinite(report.result.frequencies_cm))
    assert report.result.ped.values.shape == (definition.rank, definition.rank)
    assert (tmp_path / "csv" / "gic_gf_frequencies.csv").is_file()
    assert "ped.csv" in written
