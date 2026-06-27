from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from matrix_chem import MolecularGeometry
from matrix_chem.geometry_io import GeometryParseError, normalize_atom_symbol
from matrix_chem.zmatrix import parse_zmatrix_text, zmatrix_to_geometry


SCF_RE = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+]?\d+\.\d+)")
NORMAL_TERMINATION = "Normal termination of Gaussian"
STANDARD_ORIENTATION = "Standard orientation:"
INPUT_ORIENTATION = "Input orientation:"
POINT_GROUP_RE = re.compile(r"Full point group\s+(\S+)", flags=re.IGNORECASE)
RANK_RE = re.compile(r"NTRed=\s*(\d+).*?NRank=\s*(\d+)", flags=re.IGNORECASE)
NIMAG_RE = re.compile(r"N\s*Imag\s*=\s*(\d+)", flags=re.IGNORECASE)
READALLGIC_RE = re.compile(r"\breadallgic\b", flags=re.IGNORECASE)
PARAMETER_LINE_RE = re.compile(r"^\s*!\s*(?P<label>\S+)\s+(?P<body>.*?)\s*!\s*$")


@dataclass(frozen=True)
class GaussianLogSummary:
    path: Path
    normal_termination: bool
    scf_energies_hartree: tuple[float, ...]
    standard_orientation_count: int
    input_orientation_count: int
    scan_marker_count: int
    puckering_marker_count: int
    frequencies_cm: tuple[float, ...] = ()
    last_orientation: MolecularGeometry | None = None


@dataclass(frozen=True)
class GaussianReadAllGICLogCheck:
    path: Path
    normal_termination_count: int
    route_has_readallgic: bool
    point_groups: tuple[str, ...]
    ranks: tuple[tuple[int, int], ...]
    labels: tuple[str, ...]
    active_labels: tuple[str, ...]
    frozen_labels: tuple[str, ...]
    frequency_count: int
    imaginary_frequency_count: int
    n_imag: int | None
    optimization_completed: bool
    stationary_point: bool
    errors: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not self.errors


@dataclass(frozen=True)
class _GaussianInputBlock:
    route_lines: tuple[str, ...]
    title: str
    charge: int
    multiplicity: int
    geometry_lines: tuple[str, ...]
    tail_lines: tuple[str, ...]


def read_gaussian_input(path: Path) -> MolecularGeometry:
    """Read a Gaussian input file with Cartesian or Z-matrix geometry."""
    target = Path(path)
    block = _read_gaussian_input_block(target)
    if _geometry_looks_cartesian(block.geometry_lines):
        return _geometry_from_cartesian_block(target, block)
    return _geometry_from_zmatrix_block(target, block)


def read_gaussian_cartesian_input(path: Path) -> MolecularGeometry:
    target = Path(path)
    block = _read_gaussian_input_block(target)
    if not _geometry_looks_cartesian(block.geometry_lines) and _route_requests_zmatrix(
        block.route_lines
    ):
        return read_gaussian_zmatrix_input(target)
    return _geometry_from_cartesian_block(target, block)


def read_gaussian_zmatrix_input(path: Path) -> MolecularGeometry:
    target = Path(path)
    block = _read_gaussian_input_block(target)
    if _geometry_looks_cartesian(block.geometry_lines):
        raise GeometryParseError("Gaussian input contains Cartesian coordinates, not a Z-matrix")
    return _geometry_from_zmatrix_block(target, block)


def _read_gaussian_input_block(path: Path) -> _GaussianInputBlock:
    target = Path(path)
    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
    route_end = _route_end(lines)
    if route_end is None:
        raise GeometryParseError("Gaussian input needs a route section starting with #")
    route_lines = _route_lines(lines)
    idx = _next_nonblank(lines, route_end)
    title_start = idx
    while idx < len(lines) and lines[idx].strip():
        idx += 1
    title = " ".join(line.strip() for line in lines[title_start:idx] if line.strip())
    idx = _next_nonblank(lines, idx)
    if idx >= len(lines) or not _is_charge_multiplicity(lines[idx]):
        raise GeometryParseError("Gaussian input needs charge and multiplicity before coordinates")
    charge, multiplicity = (int(value) for value in lines[idx].split()[:2])
    idx = _next_nonblank(lines, idx + 1)
    geometry_start = idx
    while idx < len(lines) and lines[idx].strip():
        idx += 1
    geometry_lines = tuple(lines[geometry_start:idx])
    if not geometry_lines:
        raise GeometryParseError("Gaussian input contains no geometry block")
    tail_lines = tuple(lines[idx + 1 :]) if idx < len(lines) else ()
    return _GaussianInputBlock(
        route_lines=route_lines,
        title=title or target.stem,
        charge=charge,
        multiplicity=multiplicity,
        geometry_lines=geometry_lines,
        tail_lines=tail_lines,
    )


def _geometry_from_cartesian_block(path: Path, block: _GaussianInputBlock) -> MolecularGeometry:
    atoms: list[str] = []
    coords: list[list[float]] = []
    for line in block.geometry_lines:
        atom, xyz = _parse_cartesian_line(line)
        atoms.append(atom)
        coords.append(xyz)
    if not atoms:
        raise GeometryParseError("Gaussian input contains no Cartesian coordinate block")
    return MolecularGeometry(
        atoms=tuple(atoms),
        coordinates_angstrom=np.asarray(coords, dtype=float),
        comment=block.title or path.stem,
        source_format="gaussian_cartesian_input",
        source_path=path,
        charge=block.charge,
        multiplicity=block.multiplicity,
        fixed_parameters=_modredundant_fixed_patterns(block.tail_lines),
        metadata={"route": block.route_lines},
    )


def _geometry_from_zmatrix_block(path: Path, block: _GaussianInputBlock) -> MolecularGeometry:
    body_lines = [
        f"{block.charge} {block.multiplicity}",
        *block.geometry_lines,
    ]
    variable_lines = _leading_zmatrix_variable_lines(block.tail_lines)
    if variable_lines:
        body_lines.append("")
        body_lines.extend(variable_lines)
    zmat = parse_zmatrix_text("\n".join(body_lines), title=block.title)
    return zmatrix_to_geometry(zmat, source_path=path, source_format="gaussian_zmatrix_input")


def summarize_gaussian_log(path: Path) -> GaussianLogSummary:
    target = Path(path)
    text = target.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    scf = tuple(float(match.group(1)) for match in SCF_RE.finditer(text))
    orientation = _parse_last_orientation(lines, source_path=target)
    return GaussianLogSummary(
        path=target,
        normal_termination=NORMAL_TERMINATION in text,
        scf_energies_hartree=scf,
        standard_orientation_count=text.count(STANDARD_ORIENTATION),
        input_orientation_count=text.count(INPUT_ORIENTATION),
        scan_marker_count=text.lower().count("scan"),
        puckering_marker_count=text.count("QPck") + text.count("PhiP") + text.count("RPck"),
        frequencies_cm=tuple(_parse_frequencies(text)),
        last_orientation=orientation,
    )


def check_gaussian_readallgic_log(
    path: Path,
    *,
    expected_point_group: str | None = None,
    expected_rank: int | None = None,
    expected_frozen_labels: tuple[str, ...] = (),
    expected_active_labels: tuple[str, ...] = (),
    require_frequency: bool = False,
    require_no_imaginary: bool = False,
) -> GaussianReadAllGICLogCheck:
    """Validate the Gaussian-side result of a MATRIX ReadAllGIC input."""
    target = Path(path)
    text = target.read_text(encoding="utf-8", errors="replace")
    summary = summarize_gaussian_log(target)
    normal_termination_count = text.count(NORMAL_TERMINATION)
    labels, frozen_labels = _parse_parameter_labels(text)
    frozen = set(frozen_labels)
    active_labels = tuple(label for label in labels if label not in frozen)
    point_groups = tuple(match.group(1).upper() for match in POINT_GROUP_RE.finditer(text))
    ranks = tuple((int(match.group(1)), int(match.group(2))) for match in RANK_RE.finditer(text))
    n_imag = _parse_last_nimag(text)
    imaginary_frequency_count = sum(1 for frequency in summary.frequencies_cm if frequency < 0.0)

    errors: list[str] = []
    if normal_termination_count == 0:
        errors.append("Gaussian log has no normal termination")
    if not READALLGIC_RE.search(text):
        errors.append("Gaussian log does not contain readallgic")
    if expected_point_group is not None:
        expected = expected_point_group.upper()
        if expected not in point_groups:
            errors.append(f"expected point group {expected}, found {point_groups or 'none'}")
    if expected_rank is not None and (expected_rank, expected_rank) not in ranks:
        errors.append(f"expected NTRed/NRank {expected_rank}, found {ranks or 'none'}")
    missing_frozen = tuple(label for label in expected_frozen_labels if label not in frozen)
    if missing_frozen:
        errors.append(f"missing frozen labels: {', '.join(missing_frozen)}")
    missing_active = tuple(label for label in expected_active_labels if label not in labels)
    if missing_active:
        errors.append(f"missing active labels: {', '.join(missing_active)}")
    frozen_active = tuple(label for label in expected_active_labels if label in frozen)
    if frozen_active:
        errors.append(f"expected active labels are frozen: {', '.join(frozen_active)}")
    if require_frequency and not summary.frequencies_cm:
        errors.append("Gaussian log contains no vibrational frequencies")
    if require_no_imaginary:
        if n_imag is not None and n_imag != 0:
            errors.append(f"expected NImag=0, found NImag={n_imag}")
        if n_imag is None and imaginary_frequency_count:
            errors.append(f"expected no imaginary frequencies, found {imaginary_frequency_count}")

    return GaussianReadAllGICLogCheck(
        path=target,
        normal_termination_count=normal_termination_count,
        route_has_readallgic=READALLGIC_RE.search(text) is not None,
        point_groups=point_groups,
        ranks=ranks,
        labels=labels,
        active_labels=active_labels,
        frozen_labels=frozen_labels,
        frequency_count=len(summary.frequencies_cm),
        imaginary_frequency_count=imaginary_frequency_count,
        n_imag=n_imag,
        optimization_completed="Optimization completed" in text,
        stationary_point="Stationary point found" in text,
        errors=tuple(errors),
    )


def read_gaussian_log_geometry(path: Path) -> MolecularGeometry:
    """Read the last Gaussian orientation as the shared MATRIX geometry."""
    target = Path(path)
    geometry = summarize_gaussian_log(target).last_orientation
    if geometry is None:
        raise GeometryParseError("Gaussian log contains no readable orientation block")
    return geometry


def _route_end(lines: list[str]) -> int | None:
    idx = 0
    while idx < len(lines):
        if lines[idx].strip().startswith("#"):
            while idx < len(lines) and lines[idx].strip():
                idx += 1
            return idx
        idx += 1
    return None


def _route_lines(lines: list[str]) -> tuple[str, ...]:
    idx = 0
    while idx < len(lines):
        if lines[idx].strip().startswith("#"):
            route: list[str] = []
            while idx < len(lines) and lines[idx].strip():
                route.append(lines[idx].strip())
                idx += 1
            return tuple(route)
        idx += 1
    return ()


def _route_requests_zmatrix(route_lines: tuple[str, ...]) -> bool:
    route = " ".join(route_lines).lower()
    return "zmat" in route or "z-matrix" in route


def _geometry_looks_cartesian(lines: tuple[str, ...]) -> bool:
    if not lines:
        return False
    for line in lines:
        parts = line.split()
        if len(parts) < 4:
            return False
        try:
            normalize_atom_symbol(parts[0])
            _float_token(parts[1])
            _float_token(parts[2])
            _float_token(parts[3])
        except (GeometryParseError, ValueError):
            return False
    return True


def _leading_zmatrix_variable_lines(lines: tuple[str, ...]) -> tuple[str, ...]:
    variables: list[str] = []
    for raw in lines:
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped.startswith(("#", "%")):
            break
        if _is_zmatrix_variable_header(stripped) or _is_zmatrix_variable_line(stripped):
            variables.append(raw)
            continue
        break
    return tuple(variables)


def _is_zmatrix_variable_header(line: str) -> bool:
    return line.lower().rstrip(":") in {
        "variables",
        "variable",
        "constants",
        "constant",
        "parameters",
        "parameter",
    }


def _is_zmatrix_variable_line(line: str) -> bool:
    if "=" in line:
        name, value = line.split("=", 1)
    else:
        parts = line.split()
        if len(parts) != 2:
            return False
        name, value = parts
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", name.strip()):
        return False
    try:
        _float_token(value)
    except ValueError:
        return False
    return True


def _next_nonblank(lines: list[str], idx: int) -> int:
    while idx < len(lines) and not lines[idx].strip():
        idx += 1
    return idx


def _is_charge_multiplicity(line: str) -> bool:
    parts = line.split()
    if len(parts) < 2:
        return False
    try:
        int(parts[0])
        int(parts[1])
    except ValueError:
        return False
    return True


def _parse_cartesian_line(line: str) -> tuple[str, list[float]]:
    parts = line.split()
    if len(parts) < 4:
        raise GeometryParseError(f"invalid Gaussian Cartesian line: {line}")
    try:
        coords = [_float_token(parts[1]), _float_token(parts[2]), _float_token(parts[3])]
    except ValueError as exc:
        raise GeometryParseError(f"invalid Gaussian Cartesian coordinates: {line}") from exc
    return normalize_atom_symbol(parts[0]), coords


def _float_token(token: str) -> float:
    return float(token.strip().replace("D", "E").replace("d", "e"))


def _modredundant_fixed_patterns(lines: list[str] | tuple[str, ...]) -> tuple[str, ...]:
    patterns: list[str] = []
    seen: set[str] = set()
    for raw in lines:
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        parts = line.replace(",", " ").split()
        if not parts:
            continue
        kind = parts[0].upper()
        expected = {"B": 2, "A": 3, "D": 4, "O": 4, "L": 3}.get(kind)
        if expected is None:
            continue
        atoms: list[int] = []
        actions: list[str] = []
        for token in parts[1:]:
            try:
                atoms.append(int(token))
            except ValueError:
                actions.append(token.upper())
        if len(atoms) < expected or "F" not in actions:
            continue
        pattern = f"{kind}({','.join(str(atom) for atom in atoms[:expected])})"
        if pattern not in seen:
            patterns.append(pattern)
            seen.add(pattern)
    return tuple(patterns)


def _parse_frequencies(text: str) -> list[float]:
    values: list[float] = []
    for line in text.splitlines():
        if "Frequencies --" not in line:
            continue
        for token in line.split("--", 1)[1].split():
            try:
                values.append(float(token))
            except ValueError:
                continue
    return values


def _parse_parameter_labels(text: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    labels: list[str] = []
    frozen_labels: list[str] = []
    seen: set[str] = set()
    seen_frozen: set[str] = set()
    for line in text.splitlines():
        match = PARAMETER_LINE_RE.match(line)
        if not match:
            continue
        label = match.group("label")
        if not _looks_like_parameter_label(label):
            continue
        body = match.group("body")
        if label not in seen:
            labels.append(label)
            seen.add(label)
        if "frozen" in body.lower() and label not in seen_frozen:
            frozen_labels.append(label)
            seen_frozen.add(label)
    return tuple(labels), tuple(frozen_labels)


def _looks_like_parameter_label(label: str) -> bool:
    return bool(re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", label) and re.search(r"\d", label))


def _parse_last_nimag(text: str) -> int | None:
    matches = tuple(NIMAG_RE.finditer(text))
    if not matches:
        return None
    return int(matches[-1].group(1))


def _parse_last_orientation(lines: list[str], *, source_path: Path) -> MolecularGeometry | None:
    start = -1
    orientation_kind = ""
    for idx, line in enumerate(lines):
        if STANDARD_ORIENTATION in line or INPUT_ORIENTATION in line:
            start = idx
            orientation_kind = line.strip().rstrip(":")
    if start < 0:
        return None
    atoms: list[str] = []
    coords: list[list[float]] = []
    dash_count = 0
    for line in lines[start + 1 :]:
        if set(line.strip()) == {"-"}:
            dash_count += 1
            if dash_count >= 3:
                break
            continue
        if dash_count < 2:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            atom = normalize_atom_symbol(parts[1])
            xyz = [float(parts[3]), float(parts[4]), float(parts[5])]
        except ValueError:
            continue
        atoms.append(atom)
        coords.append(xyz)
    if not atoms:
        return None
    return MolecularGeometry(
        atoms=tuple(atoms),
        coordinates_angstrom=np.asarray(coords, dtype=float),
        comment=f"Gaussian {orientation_kind}",
        source_format="gaussian_log_orientation",
        source_path=source_path,
        metadata={"orientation": orientation_kind},
    )
