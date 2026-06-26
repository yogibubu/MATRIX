from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from oracle_core import (
    key_value_section_lines,
    parse_key_value_section,
    read_sectioned_lines,
    replace_section,
    section_content,
)


ORACLE_XYZ_CARTESIAN_HESSIAN_SCHEMA = "oracle.xyz.cartesian_hessian.v1"
ORACLE_XYZ_NORMAL_MODES_SCHEMA = "oracle.xyz.normal_modes.v1"
ORACLE_XYZ_QFF_SCHEMA = "oracle.xyz.qff.v1"


@dataclass(frozen=True)
class CartesianHessianSection:
    atomic_numbers: tuple[int, ...]
    cartesian_coordinates_bohr: np.ndarray
    masses_amu: tuple[float, ...]
    cartesian_hessian: np.ndarray
    harmonic_frequencies_cm: tuple[float, ...] = ()
    source: str = ""
    schema: str = ORACLE_XYZ_CARTESIAN_HESSIAN_SCHEMA

    def __post_init__(self) -> None:
        natoms = len(self.atomic_numbers)
        coords = np.asarray(self.cartesian_coordinates_bohr, dtype=float)
        hessian = np.asarray(self.cartesian_hessian, dtype=float)
        if coords.shape != (natoms, 3):
            raise ValueError(f"coordinates must have shape ({natoms}, 3)")
        if len(self.masses_amu) != natoms:
            raise ValueError(f"masses must have length {natoms}")
        expected = (3 * natoms, 3 * natoms)
        if hessian.shape != expected:
            raise ValueError(f"Cartesian Hessian must have shape {expected}")
        if not np.allclose(hessian, hessian.T):
            raise ValueError("Cartesian Hessian must be symmetric")
        object.__setattr__(self, "atomic_numbers", tuple(int(value) for value in self.atomic_numbers))
        object.__setattr__(self, "cartesian_coordinates_bohr", coords)
        object.__setattr__(self, "masses_amu", tuple(float(value) for value in self.masses_amu))
        object.__setattr__(self, "cartesian_hessian", hessian)
        object.__setattr__(
            self,
            "harmonic_frequencies_cm",
            tuple(float(value) for value in self.harmonic_frequencies_cm),
        )


@dataclass(frozen=True)
class NormalModesSection:
    frequencies_cm: tuple[float, ...]
    modes: np.ndarray
    source: str = ""
    schema: str = ORACLE_XYZ_NORMAL_MODES_SCHEMA

    def __post_init__(self) -> None:
        modes = np.asarray(self.modes, dtype=float)
        if modes.ndim != 2:
            raise ValueError("normal modes must be a mode x Cartesian-coordinate matrix")
        if self.frequencies_cm and modes.shape[0] != len(self.frequencies_cm):
            raise ValueError("normal-mode row count must match frequency count")
        object.__setattr__(self, "frequencies_cm", tuple(float(value) for value in self.frequencies_cm))
        object.__setattr__(self, "modes", modes)


@dataclass(frozen=True)
class QFFSection:
    harmonic_frequencies_cm: tuple[float, ...]
    anharmonic_frequencies_cm: tuple[float, ...] = ()
    cubic_cm: dict[tuple[int, int, int], float] | None = None
    quartic_cm: dict[tuple[int, int, int, int], float] | None = None
    source: str = ""
    schema: str = ORACLE_XYZ_QFF_SCHEMA

    def __post_init__(self) -> None:
        cubic = {
            tuple(sorted(int(item) for item in key)): float(value)
            for key, value in (self.cubic_cm or {}).items()
        }
        quartic = {
            tuple(sorted(int(item) for item in key)): float(value)
            for key, value in (self.quartic_cm or {}).items()
        }
        object.__setattr__(
            self,
            "harmonic_frequencies_cm",
            tuple(float(value) for value in self.harmonic_frequencies_cm),
        )
        object.__setattr__(
            self,
            "anharmonic_frequencies_cm",
            tuple(float(value) for value in self.anharmonic_frequencies_cm),
        )
        object.__setattr__(self, "cubic_cm", cubic)
        object.__setattr__(self, "quartic_cm", quartic)


def cartesian_hessian_section_from_hessian_input(input_data, *, source: str | None = None) -> CartesianHessianSection:
    return CartesianHessianSection(
        atomic_numbers=tuple(int(value) for value in input_data.atomic_numbers),
        cartesian_coordinates_bohr=np.asarray(input_data.cartesian_coordinates_bohr, dtype=float),
        masses_amu=tuple(float(value) for value in input_data.masses_amu),
        cartesian_hessian=np.asarray(input_data.cartesian_hessian, dtype=float),
        harmonic_frequencies_cm=tuple(float(value) for value in input_data.harmonic_frequencies_cm),
        source=input_data.source if source is None else source,
    )


def hessian_input_from_cartesian_hessian_section(section: CartesianHessianSection):
    from oracle_gf import HessianInput

    data = HessianInput(
        atomic_numbers=np.asarray(section.atomic_numbers, dtype=int),
        cartesian_coordinates_bohr=np.asarray(section.cartesian_coordinates_bohr, dtype=float),
        masses_amu=np.asarray(section.masses_amu, dtype=float),
        cartesian_hessian=np.asarray(section.cartesian_hessian, dtype=float),
        harmonic_frequencies_cm=np.asarray(section.harmonic_frequencies_cm, dtype=float),
        source=section.source or "xyzin-cartesian-hessian",
    )
    data.validate()
    return data


def hessian_input_from_xyzin(path: Path | str):
    return hessian_input_from_cartesian_hessian_section(read_cartesian_hessian_section(Path(path)))


def normal_modes_section_from_arrays(
    frequencies_cm,
    modes,
    *,
    source: str = "",
    coordinate_count: int | None = None,
) -> NormalModesSection:
    values = np.asarray(modes, dtype=float).reshape(-1)
    freqs = tuple(float(value) for value in np.asarray(frequencies_cm, dtype=float).reshape(-1))
    if coordinate_count is None:
        coordinate_count = len(values) // len(freqs) if freqs else len(values)
    if coordinate_count <= 0:
        raise ValueError("normal-mode coordinate count must be positive")
    if values.size % coordinate_count != 0:
        raise ValueError("normal-mode flat array is not divisible by coordinate count")
    matrix = values.reshape((values.size // coordinate_count, coordinate_count))
    if freqs and matrix.shape[0] != len(freqs):
        matrix = matrix[: len(freqs), :]
    return NormalModesSection(frequencies_cm=freqs, modes=matrix, source=source)


def qff_section_from_anharmonic_input(input_data, *, source: str | None = None) -> QFFSection:
    return QFFSection(
        harmonic_frequencies_cm=tuple(float(value) for value in input_data.harmonic_frequencies_cm),
        anharmonic_frequencies_cm=tuple(float(value) for value in input_data.anharmonic_frequencies_cm),
        cubic_cm=dict(input_data.cubic_cm),
        quartic_cm=dict(input_data.quartic_cm),
        source=input_data.source if source is None else source,
    )


def anharmonic_input_from_qff_section(section: QFFSection):
    from oracle_vpt2_vci import AnharmonicInput

    data = AnharmonicInput(
        harmonic_frequencies_cm=np.asarray(section.harmonic_frequencies_cm, dtype=float),
        anharmonic_frequencies_cm=np.asarray(section.anharmonic_frequencies_cm, dtype=float),
        cubic_cm=dict(section.cubic_cm or {}),
        quartic_cm=dict(section.quartic_cm or {}),
        source=section.source or "xyzin-qff",
    )
    data.validate()
    return data


def qff_section_from_quartic_force_field(force_field, *, source: str = "") -> QFFSection:
    return QFFSection(
        harmonic_frequencies_cm=tuple(float(value) for value in force_field.harmonic_frequencies_cm),
        cubic_cm=dict(force_field.cubic_cm),
        quartic_cm=dict(force_field.quartic_cm),
        source=source,
    )


def quartic_force_field_from_qff_section(section: QFFSection):
    from oracle_vpt2_vci import QuarticForceField

    frequencies = (
        section.anharmonic_frequencies_cm
        if section.anharmonic_frequencies_cm
        else section.harmonic_frequencies_cm
    )
    return QuarticForceField(
        np.asarray(frequencies, dtype=float),
        dict(section.cubic_cm or {}),
        dict(section.quartic_cm or {}),
    )


def cartesian_hessian_section_lines(section: CartesianHessianSection) -> list[str]:
    lines = key_value_section_lines(
        ORACLE_XYZ_CARTESIAN_HESSIAN_SCHEMA,
        {
            "SOURCE": section.source or None,
            "NATOMS": len(section.atomic_numbers),
            "ATOMIC_NUMBERS": " ".join(str(value) for value in section.atomic_numbers),
            "MASSES_AMU": _format_float_list(section.masses_amu),
            "HARMONIC_FREQ_CM1": _format_float_list(section.harmonic_frequencies_cm),
        },
        key_order=("SOURCE", "NATOMS", "ATOMIC_NUMBERS", "MASSES_AMU", "HARMONIC_FREQ_CM1"),
    )
    lines.append("[COORDINATES_BOHR]")
    for idx, (x, y, z) in enumerate(section.cartesian_coordinates_bohr, start=1):
        lines.append(f"{idx:d} {_format_float(x)} {_format_float(y)} {_format_float(z)}")
    lines.append("[HESSIAN_LOWER_AU]")
    for idx, value in enumerate(_lower_from_symmetric(section.cartesian_hessian), start=1):
        lines.append(f"{idx:d} {_format_float(value)}")
    return lines


def parse_cartesian_hessian_section(lines: list[str] | tuple[str, ...]) -> CartesianHessianSection:
    raw_lines = list(lines)
    values = parse_key_value_section(raw_lines)
    atomic_numbers = tuple(int(value) for value in _number_list(values.get("ATOMIC_NUMBERS")))
    natoms = int(float(values.get("NATOMS", len(atomic_numbers))))
    if len(atomic_numbers) != natoms:
        raise ValueError("CARTESIAN_HESSIAN atomic number count does not match NATOMS")
    coords = _read_coordinate_rows(_section_block(raw_lines, "COORDINATES_BOHR"), natoms=natoms)
    lower = _read_indexed_values(_section_block(raw_lines, "HESSIAN_LOWER_AU"))
    return CartesianHessianSection(
        atomic_numbers=atomic_numbers,
        cartesian_coordinates_bohr=coords,
        masses_amu=tuple(_number_list(values.get("MASSES_AMU"))),
        cartesian_hessian=_symmetric_from_lower(np.asarray(lower, dtype=float)),
        harmonic_frequencies_cm=tuple(_number_list(values.get("HARMONIC_FREQ_CM1"))),
        source=values.get("SOURCE", ""),
        schema=values.get("SCHEMA", ORACLE_XYZ_CARTESIAN_HESSIAN_SCHEMA),
    )


def write_cartesian_hessian_section(path: Path | str, section: CartesianHessianSection) -> None:
    replace_section(Path(path), "CARTESIAN_HESSIAN", cartesian_hessian_section_lines(section))


def read_cartesian_hessian_section(path: Path | str) -> CartesianHessianSection:
    content = section_content(read_sectioned_lines(Path(path)), "CARTESIAN_HESSIAN")
    if not content:
        raise ValueError("missing #CARTESIAN_HESSIAN section")
    return parse_cartesian_hessian_section(content)


def normal_modes_section_lines(section: NormalModesSection) -> list[str]:
    lines = key_value_section_lines(
        ORACLE_XYZ_NORMAL_MODES_SCHEMA,
        {
            "SOURCE": section.source or None,
            "MODE_COUNT": section.modes.shape[0],
            "COORDINATE_COUNT": section.modes.shape[1],
            "FREQ_CM1": _format_float_list(section.frequencies_cm),
        },
        key_order=("SOURCE", "MODE_COUNT", "COORDINATE_COUNT", "FREQ_CM1"),
    )
    lines.append("[MODES]")
    for idx, value in enumerate(section.modes.reshape(-1), start=1):
        lines.append(f"{idx:d} {_format_float(value)}")
    return lines


def parse_normal_modes_section(lines: list[str] | tuple[str, ...]) -> NormalModesSection:
    raw_lines = list(lines)
    values = parse_key_value_section(raw_lines)
    mode_count = int(float(values.get("MODE_COUNT", "0")))
    coordinate_count = int(float(values.get("COORDINATE_COUNT", "0")))
    flat = np.asarray(_read_indexed_values(_section_block(raw_lines, "MODES")), dtype=float)
    if mode_count <= 0 or coordinate_count <= 0:
        raise ValueError("NORMAL_MODES needs MODE_COUNT and COORDINATE_COUNT")
    if flat.size != mode_count * coordinate_count:
        raise ValueError("NORMAL_MODES flat value count does not match declared shape")
    return NormalModesSection(
        frequencies_cm=tuple(_number_list(values.get("FREQ_CM1"))),
        modes=flat.reshape((mode_count, coordinate_count)),
        source=values.get("SOURCE", ""),
        schema=values.get("SCHEMA", ORACLE_XYZ_NORMAL_MODES_SCHEMA),
    )


def write_normal_modes_section(path: Path | str, section: NormalModesSection) -> None:
    replace_section(Path(path), "NORMAL_MODES", normal_modes_section_lines(section))


def read_normal_modes_section(path: Path | str) -> NormalModesSection:
    content = section_content(read_sectioned_lines(Path(path)), "NORMAL_MODES")
    if not content:
        raise ValueError("missing #NORMAL_MODES section")
    return parse_normal_modes_section(content)


def qff_section_lines(section: QFFSection) -> list[str]:
    lines = key_value_section_lines(
        ORACLE_XYZ_QFF_SCHEMA,
        {
            "SOURCE": section.source or None,
            "MODE_COUNT": len(section.harmonic_frequencies_cm),
            "HARMONIC_FREQ_CM1": _format_float_list(section.harmonic_frequencies_cm),
            "ANHARMONIC_FREQ_CM1": _format_float_list(section.anharmonic_frequencies_cm),
        },
        key_order=("SOURCE", "MODE_COUNT", "HARMONIC_FREQ_CM1", "ANHARMONIC_FREQ_CM1"),
    )
    if section.cubic_cm:
        lines.append("[CUBIC_CM]")
        for key, value in sorted(section.cubic_cm.items()):
            i, j, k = (item + 1 for item in key)
            lines.append(f"{i:d} {j:d} {k:d} {_format_float(value)}")
    if section.quartic_cm:
        lines.append("[QUARTIC_CM]")
        for key, value in sorted(section.quartic_cm.items()):
            i, j, k, l = (item + 1 for item in key)
            lines.append(f"{i:d} {j:d} {k:d} {l:d} {_format_float(value)}")
    return lines


def parse_qff_section(lines: list[str] | tuple[str, ...]) -> QFFSection:
    raw_lines = list(lines)
    values = parse_key_value_section(raw_lines)
    return QFFSection(
        harmonic_frequencies_cm=tuple(_number_list(values.get("HARMONIC_FREQ_CM1"))),
        anharmonic_frequencies_cm=tuple(_number_list(values.get("ANHARMONIC_FREQ_CM1"))),
        cubic_cm=_read_force_terms(_section_block(raw_lines, "CUBIC_CM"), order=3),
        quartic_cm=_read_force_terms(_section_block(raw_lines, "QUARTIC_CM"), order=4),
        source=values.get("SOURCE", ""),
        schema=values.get("SCHEMA", ORACLE_XYZ_QFF_SCHEMA),
    )


def write_qff_section(path: Path | str, section: QFFSection) -> None:
    replace_section(Path(path), "QFF", qff_section_lines(section))


def read_qff_section(path: Path | str) -> QFFSection:
    content = section_content(read_sectioned_lines(Path(path)), "QFF")
    if not content:
        raise ValueError("missing #QFF section")
    return parse_qff_section(content)


def _section_block(lines: list[str], name: str) -> list[str]:
    header = f"[{name.upper()}]"
    out: list[str] = []
    active = False
    for raw in lines:
        text = raw.strip()
        if text.upper() == header:
            active = True
            continue
        if active and text.startswith("[") and text.endswith("]"):
            break
        if active and text:
            out.append(text)
    return out


def _read_coordinate_rows(lines: list[str], *, natoms: int) -> np.ndarray:
    rows: list[list[float]] = []
    for line in lines:
        values = _number_list(line)
        if len(values) >= 4:
            rows.append([float(values[-3]), float(values[-2]), float(values[-1])])
    coords = np.asarray(rows, dtype=float)
    if coords.shape != (natoms, 3):
        raise ValueError("COORDINATES_BOHR row count does not match NATOMS")
    return coords


def _read_indexed_values(lines: list[str]) -> list[float]:
    values: list[float] = []
    for line in lines:
        numbers = _number_list(line)
        if numbers:
            values.append(float(numbers[-1]))
    return values


def _read_force_terms(lines: list[str], *, order: int) -> dict[tuple[int, ...], float]:
    out: dict[tuple[int, ...], float] = {}
    for line in lines:
        numbers = _number_list(line)
        if len(numbers) < order + 1:
            continue
        key = tuple(sorted(int(value) - 1 for value in numbers[:order]))
        out[key] = float(numbers[order])
    return out


def _number_list(text: str | None) -> list[float]:
    if not text:
        return []
    values: list[float] = []
    for token in str(text).replace(",", " ").split():
        try:
            values.append(float(token.replace("D", "E").replace("d", "e")))
        except ValueError:
            continue
    return values


def _format_float(value: float) -> str:
    return f"{float(value):.12g}"


def _format_float_list(values) -> str | None:
    values_tuple = tuple(values)
    if not values_tuple:
        return None
    return " ".join(_format_float(float(value)) for value in values_tuple)


def _lower_from_symmetric(matrix: np.ndarray) -> np.ndarray:
    mat = np.asarray(matrix, dtype=float)
    values: list[float] = []
    for i in range(mat.shape[0]):
        for j in range(i + 1):
            values.append(float(mat[i, j]))
    return np.asarray(values, dtype=float)


def _symmetric_from_lower(lower: np.ndarray) -> np.ndarray:
    n_float = (np.sqrt(8 * len(lower) + 1) - 1) / 2
    n = int(round(n_float))
    if n * (n + 1) // 2 != len(lower):
        raise ValueError("packed lower-triangular Hessian has invalid length")
    mat = np.zeros((n, n), dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i + 1):
            mat[i, j] = lower[idx]
            mat[j, i] = lower[idx]
            idx += 1
    return mat
