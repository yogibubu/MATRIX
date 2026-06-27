from pathlib import Path

from matrix_gaussian import check_gaussian_readallgic_log


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
