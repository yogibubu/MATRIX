from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
import re

from oracle_core import (
    key_value_section_lines,
    normalize_key,
    parse_key_value_section,
    read_sectioned_lines,
    replace_section,
    section_content,
)


ORACLE_XYZ_ROTATIONAL_SCHEMA = "oracle.xyz.rotational.v1"
ORACLE_XYZ_VIBRATIONAL_SCHEMA = "oracle.xyz.vibrational.v1"
ORACLE_XYZ_DELTABVIB_SCHEMA = "oracle.xyz.deltabvib.v1"
ORACLE_XYZ_CORIOLIS_SCHEMA = "oracle.xyz.coriolis.v1"
ORACLE_XYZ_QCENT_SCHEMA = "oracle.xyz.qcent.v1"

MERLINO_XYZIN_ROTATIONAL_SCHEMA = "merlino.xyzin.rotational.v1"
MERLINO_XYZIN_VIBRATIONAL_SCHEMA = "merlino.xyzin.vibrational.v1"


@dataclass(frozen=True)
class RotationalSection:
    rotor_type: str = ""
    representation: str = ""
    point_group: str = ""
    symmetry_number: int | None = None
    temperature_K: float | None = None
    pressure_atm: float | None = None
    A_MHz: float | None = None
    B_MHz: float | None = None
    C_MHz: float | None = None
    q_rot: float | None = None
    delta_vib_MHz: tuple[float | None, float | None, float | None] | None = None
    schema: str = ORACLE_XYZ_ROTATIONAL_SCHEMA


@dataclass(frozen=True)
class VibrationalSection:
    linear: bool | None = None
    nvib: int | None = None
    n_imag_like: int | None = None
    symmetry_group: str = ""
    frequencies_cm1: tuple[float, ...] = ()
    ir_intensities_km_mol: tuple[float, ...] = ()
    chi_cm1: tuple[tuple[int, int, float], ...] = ()
    schema: str = ORACLE_XYZ_VIBRATIONAL_SCHEMA


def parse_rotational_section(lines: Iterable[str]) -> RotationalSection:
    values = parse_key_value_section(lines)
    return RotationalSection(
        rotor_type=values.get("ROTOR_TYPE", ""),
        representation=values.get("REPRESENTATION", ""),
        point_group=values.get("POINT_GROUP", ""),
        symmetry_number=_optional_int(values.get("SYMM_NUMBER") or values.get("SIGMA")),
        temperature_K=_optional_float(values.get("T_K")),
        pressure_atm=_optional_float(values.get("P_ATM")),
        A_MHz=_optional_float(values.get("A_MHZ")),
        B_MHz=_optional_float(values.get("B_MHZ")),
        C_MHz=_optional_float(values.get("C_MHZ")),
        q_rot=_optional_float(values.get("Q_ROT")),
        delta_vib_MHz=_parse_delta_vib(values),
        schema=values.get("SCHEMA", ORACLE_XYZ_ROTATIONAL_SCHEMA),
    )


def rotational_section_lines(section: RotationalSection) -> list[str]:
    values: dict[str, object] = {
        "ROTOR_TYPE": section.rotor_type or None,
        "REPRESENTATION": section.representation or None,
        "POINT_GROUP": section.point_group or None,
        "SYMM_NUMBER": section.symmetry_number,
        "T_K": _format_float(section.temperature_K),
        "P_ATM": _format_float(section.pressure_atm),
        "A_MHZ": _format_float(section.A_MHz),
        "B_MHZ": _format_float(section.B_MHz),
        "C_MHZ": _format_float(section.C_MHz),
        "Q_ROT": _format_float(section.q_rot),
    }
    if section.delta_vib_MHz is not None:
        a, b, c = section.delta_vib_MHz
        values.update(
            {
                "DVIBA_MHZ": _format_float(a),
                "DVIBB_MHZ": _format_float(b),
                "DVIBC_MHZ": _format_float(c),
            }
        )
    return key_value_section_lines(
        ORACLE_XYZ_ROTATIONAL_SCHEMA,
        values,
        key_order=(
            "ROTOR_TYPE",
            "REPRESENTATION",
            "POINT_GROUP",
            "SYMM_NUMBER",
            "T_K",
            "P_ATM",
            "A_MHZ",
            "B_MHZ",
            "C_MHZ",
            "DVIBA_MHZ",
            "DVIBB_MHZ",
            "DVIBC_MHZ",
            "Q_ROT",
        ),
    )


def read_rotational_section(path: Path) -> RotationalSection:
    return parse_rotational_section(section_content(read_sectioned_lines(Path(path)), "ROTATIONAL"))


def write_rotational_section(path: Path, section: RotationalSection) -> None:
    replace_section(Path(path), "ROTATIONAL", rotational_section_lines(section))


def parse_vibrational_section(lines: Iterable[str]) -> VibrationalSection:
    raw_lines = list(lines)
    values = parse_key_value_section(raw_lines)
    return VibrationalSection(
        linear=_optional_bool(values.get("LINEAR")),
        nvib=_optional_int(values.get("NVIB")),
        n_imag_like=_optional_int(values.get("N_IMAG_LIKE")),
        symmetry_group=values.get("SYMMETRY_GROUP", ""),
        frequencies_cm1=tuple(_number_list(values.get("FREQ_CM1") or values.get("FREQUENCIES"))),
        ir_intensities_km_mol=tuple(
            _number_list(values.get("IR_INTEN_KM_MOL") or values.get("IR_INTEN"))
        ),
        chi_cm1=tuple(_parse_chi_block(raw_lines)),
        schema=values.get("SCHEMA", ORACLE_XYZ_VIBRATIONAL_SCHEMA),
    )


def vibrational_section_lines(section: VibrationalSection) -> list[str]:
    values: dict[str, object] = {
        "LINEAR": None if section.linear is None else int(bool(section.linear)),
        "NVIB": section.nvib,
        "N_IMAG_LIKE": section.n_imag_like,
        "SYMMETRY_GROUP": section.symmetry_group or None,
        "FREQ_CM1": _format_float_list(section.frequencies_cm1),
        "IR_INTEN_KM_MOL": _format_float_list(section.ir_intensities_km_mol),
    }
    lines = key_value_section_lines(
        ORACLE_XYZ_VIBRATIONAL_SCHEMA,
        values,
        key_order=(
            "LINEAR",
            "NVIB",
            "N_IMAG_LIKE",
            "SYMMETRY_GROUP",
            "FREQ_CM1",
            "IR_INTEN_KM_MOL",
        ),
    )
    if section.chi_cm1:
        lines.append("CHI_CM1 = [")
        lines.extend(f"{i:d} {j:d} {value:.8f}" for i, j, value in section.chi_cm1)
        lines.append("]")
    return lines


def read_vibrational_section(path: Path) -> VibrationalSection:
    content = section_content(read_sectioned_lines(Path(path)), "VIBRATIONAL")
    return parse_vibrational_section(content)


def write_vibrational_section(path: Path, section: VibrationalSection) -> None:
    replace_section(Path(path), "VIBRATIONAL", vibrational_section_lines(section))


def _parse_delta_vib(
    values: dict[str, str],
) -> tuple[float | None, float | None, float | None] | None:
    a = _optional_float(values.get("DVIBA_MHZ") or values.get("DVIB_A_MHZ"))
    b = _optional_float(values.get("DVIBB_MHZ") or values.get("DVIB_B_MHZ"))
    c = _optional_float(values.get("DVIBC_MHZ") or values.get("DVIB_C_MHZ"))
    if a is None and b is None and c is None:
        return None
    return a, b, c


def _parse_chi_block(lines: list[str]) -> list[tuple[int, int, float]]:
    out: list[tuple[int, int, float]] = []
    in_chi = False
    for raw in lines:
        line = raw.strip()
        upper = line.upper()
        if upper.startswith("CHI_CM1"):
            in_chi = True
            continue
        if in_chi and "]" in line:
            break
        if not in_chi or not line:
            continue
        values = _number_list(line)
        if len(values) >= 3:
            out.append((int(values[0]), int(values[1]), float(values[2])))
    return out


def _number_list(text: str | None) -> list[float]:
    if not text:
        return []
    return [
        float(item.replace("D", "E").replace("d", "e"))
        for item in re.findall(r"[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?", text)
    ]


def _optional_float(text: str | None) -> float | None:
    values = _number_list(text)
    return values[0] if values else None


def _optional_int(text: str | None) -> int | None:
    value = _optional_float(text)
    return None if value is None else int(value)


def _optional_bool(text: str | None) -> bool | None:
    if text is None:
        return None
    key = normalize_key(text)
    if key in {"1", "TRUE", "YES", "Y"}:
        return True
    if key in {"0", "FALSE", "NO", "N"}:
        return False
    value = _optional_int(text)
    return None if value is None else bool(value)


def _format_float(value: float | None) -> str | None:
    return None if value is None else f"{float(value):.8g}"


def _format_float_list(values: tuple[float, ...]) -> str | None:
    if not values:
        return None
    return " ".join(f"{float(value):.8f}" for value in values)
