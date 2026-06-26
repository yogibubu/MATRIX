from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


_HEADER_RE = re.compile(
    r"^(?P<label>.+?)\s+(?P<kind>[IRC])(?:\s+N=\s*(?P<count>\d+)|\s+(?P<scalar>[-+]?\d+))\s*$"
)


@dataclass(frozen=True)
class FCHKData:
    """Gaussian FCHK adapter output for harmonic and anharmonic consumers."""

    atomic_numbers: np.ndarray
    cartesian_coordinates_bohr: np.ndarray
    masses_amu: np.ndarray
    cartesian_hessian_lower: np.ndarray
    harmonic_frequencies_cm: np.ndarray
    anharmonic_frequencies_cm: np.ndarray
    anharmonic_e2: np.ndarray
    normal_modes: np.ndarray

    def to_hessian_input(self):
        from oracle_gf import HessianInput

        hessian = lower_to_symmetric(self.cartesian_hessian_lower)
        data = HessianInput(
            atomic_numbers=self.atomic_numbers,
            cartesian_coordinates_bohr=self.cartesian_coordinates_bohr,
            masses_amu=self.masses_amu,
            cartesian_hessian=hessian,
            harmonic_frequencies_cm=self.harmonic_frequencies_cm,
            source="gaussian-fchk",
        )
        data.validate()
        return data


def read_gaussian_fchk(path: Path) -> FCHKData:
    """Read harmonic Hessian and Gaussian anharmonic arrays from an FCHK file."""
    blocks = _read_fchk_blocks(Path(path))
    atomic_numbers = _first_array(blocks, "Atomic numbers").astype(int)
    coords = _first_array(blocks, "Current cartesian coordinates").reshape((-1, 3))
    masses = _first_array(
        blocks,
        "Real atomic weights",
        "Atomic masses",
        "Vib-AtMass",
        "Anharmonic Vib-AtMass",
    )
    hessian = _first_array(blocks, "Cartesian Force Constants")
    vib_e2 = _first_array(blocks, "Vib-E2") if "Vib-E2" in blocks else np.array((), dtype=float)
    anh_e2 = (
        _first_array(blocks, "Anharmonic Vib-E2")
        if "Anharmonic Vib-E2" in blocks
        else np.array((), dtype=float)
    )
    modes = (
        _first_array(blocks, "Anharmonic Vib-Modes", "Vib-Modes")
        if ("Anharmonic Vib-Modes" in blocks or "Vib-Modes" in blocks)
        else np.array((), dtype=float)
    )

    harmonic = (
        vib_e2[: int(blocks.get("Vib-NDim", len(masses) * 3))]
        if vib_e2.size
        else np.array((), dtype=float)
    )
    n_anh = int(blocks.get("Anharmonic Vib-NDim", harmonic.size))
    anharmonic = anh_e2[:n_anh] if anh_e2.size else np.array((), dtype=float)
    return FCHKData(
        atomic_numbers=atomic_numbers,
        cartesian_coordinates_bohr=coords,
        masses_amu=masses,
        cartesian_hessian_lower=hessian,
        harmonic_frequencies_cm=harmonic,
        anharmonic_frequencies_cm=anharmonic,
        anharmonic_e2=anh_e2,
        normal_modes=modes,
    )


def hessian_input_from_gaussian_fchk(path: Path):
    """Gaussian FCHK adapter: return the canonical ORACLE Hessian input."""
    return read_gaussian_fchk(path).to_hessian_input()


def lower_to_symmetric(lower: np.ndarray) -> np.ndarray:
    """Convert Gaussian lower-triangular packed storage to a full matrix."""
    n_float = (np.sqrt(8 * len(lower) + 1) - 1) / 2
    n = int(round(n_float))
    if n * (n + 1) // 2 != len(lower):
        raise ValueError("Packed lower-triangular array has an invalid length")
    mat = np.zeros((n, n), dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i + 1):
            mat[i, j] = lower[idx]
            mat[j, i] = lower[idx]
            idx += 1
    return mat


def _read_fchk_blocks(path: Path) -> dict[str, np.ndarray | int | float | str]:
    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    blocks: dict[str, np.ndarray | int | float | str] = {}
    i = 0
    while i < len(lines):
        match = _HEADER_RE.match(lines[i])
        if match is None:
            i += 1
            continue
        label = match.group("label").strip()
        kind = match.group("kind")
        count = match.group("count")
        scalar = match.group("scalar")
        if count is None:
            if scalar is not None:
                blocks[label] = int(scalar) if kind == "I" else float(scalar)
            i += 1
            continue

        nvalues = int(count)
        raw: list[str] = []
        i += 1
        while len(raw) < nvalues and i < len(lines):
            raw.extend(lines[i].split())
            i += 1
        if kind == "I":
            blocks[label] = np.array([int(x) for x in raw[:nvalues]], dtype=int)
        elif kind == "R":
            blocks[label] = np.array([float(x.replace("D", "E")) for x in raw[:nvalues]], dtype=float)
        else:
            blocks[label] = " ".join(raw[:nvalues])
    return blocks


def _first_array(blocks: dict[str, object], *labels: str) -> np.ndarray:
    for label in labels:
        value = blocks.get(label)
        if isinstance(value, np.ndarray):
            return value.astype(float, copy=False)
    raise ValueError(f"FCHK block not found: {' / '.join(labels)}")
