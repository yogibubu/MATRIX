from __future__ import annotations

from oracle_core import section_content
from oracle_rovib import (
    ORACLE_XYZ_ROTATIONAL_SCHEMA,
    ORACLE_XYZ_VIBRATIONAL_SCHEMA,
    RotationalSection,
    VibrationalSection,
    parse_rotational_section,
    parse_vibrational_section,
    read_rotational_section,
    read_vibrational_section,
    rovib_summary_lines,
    summarize_xyzin,
    write_rotational_section,
    write_vibrational_section,
)
from oracle_thermo import ORACLE_XYZ_THERMO_SCHEMA


def test_rovib_and_thermo_schema_constants_are_oracle_xyz_sections():
    assert ORACLE_XYZ_ROTATIONAL_SCHEMA == "oracle.xyz.rotational.v1"
    assert ORACLE_XYZ_VIBRATIONAL_SCHEMA == "oracle.xyz.vibrational.v1"
    assert ORACLE_XYZ_THERMO_SCHEMA == "oracle.xyz.thermo.v1"


def test_rotational_section_accepts_merlino_keys():
    section = parse_rotational_section(
        [
            "rotor_type = asymmetric_top_quasi_prolate",
            "Point Group = C2v",
            "Symm. Number = 2",
            "A_MHz = 1000.0",
            "B_MHz = 800.0",
            "C_MHz = 600.0",
            "DVibA_MHz=1.0",
            "DVibB_MHz=2.0",
            "DVibC_MHz=3.0",
            "Q_rot = 10.5",
        ]
    )

    assert section.point_group == "C2v"
    assert section.symmetry_number == 2
    assert section.A_MHz == 1000.0
    assert section.delta_vib_MHz == (1.0, 2.0, 3.0)
    assert section.q_rot == 10.5


def test_rotational_section_writer_preserves_other_sections(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text("1\nh\nH 0 0 0\n\n#GIC\nSCHEMA oracle.xyz.gic.v1\n", encoding="utf-8")

    write_rotational_section(
        path,
        RotationalSection(
            rotor_type="linear_top",
            B_MHz=42.0,
            symmetry_number=1,
            temperature_K=298.15,
            pressure_atm=1.0,
        ),
    )
    parsed = read_rotational_section(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    assert parsed.B_MHz == 42.0
    assert section_content(lines, "GIC")[0] == "SCHEMA oracle.xyz.gic.v1"


def test_vibrational_section_reads_frequencies_and_chi_block():
    section = parse_vibrational_section(
        [
            "linear = 0",
            "nvib = 3",
            "n_imag_like = 1",
            "symmetry_group = Cs",
            "freq_cm1 = 100.0 200.0 300.0",
            "ir_inten_km_mol = 1.0 2.0 3.0",
            "chi_cm1 = [",
            "1 1 -0.5",
            "2 1 0.1",
            "]",
        ]
    )

    assert section.linear is False
    assert section.nvib == 3
    assert section.frequencies_cm1 == (100.0, 200.0, 300.0)
    assert section.ir_intensities_km_mol == (1.0, 2.0, 3.0)
    assert section.chi_cm1 == ((1, 1, -0.5), (2, 1, 0.1))


def test_vibrational_section_writer_round_trips(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text("1\nh\nH 0 0 0\n", encoding="utf-8")

    write_vibrational_section(
        path,
        VibrationalSection(
            linear=False,
            nvib=2,
            frequencies_cm1=(100.0, 200.0),
            chi_cm1=((1, 1, -0.5),),
        ),
    )
    parsed = read_vibrational_section(path)

    assert parsed.nvib == 2
    assert parsed.frequencies_cm1 == (100.0, 200.0)
    assert parsed.chi_cm1 == ((1, 1, -0.5),)


def test_rovib_summary_reads_standalone_xyzin(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#BASIC",
                "CHARGE 0",
                "SPIN_MULTIPLICITY 1",
                "POINT_GROUP C1",
                "",
                "#ROTATIONAL",
                "A_MHz = 1000.0",
                "B_MHz = 900.0",
                "C_MHz = 800.0",
                "",
                "#VIBRATIONAL",
                "freq_cm1 = 100.0 200.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = summarize_xyzin(path)
    lines = rovib_summary_lines(summary)

    assert summary.basic.point_group == "C1"
    assert summary.rotational is not None
    assert summary.vibrational is not None
    assert "rotational: A=1000MHz B=900MHz C=800MHz" in lines
    assert "vibrational: 2 frequencies" in lines
