from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile

import numpy as np

from oracle_chem import preprocess_to_enriched_xyz, write_validation_section
from oracle_chem.topology.elements import atomic_symbol
from oracle_gicforge import (
    GICDefinition,
    build_gic_b_matrix,
    read_gic_definition_from_xyzin,
    write_gicforge_build_sections,
)

from .harmonic import solve_wilson_gf
from .models import HessianInput


BOHR_TO_ANGSTROM = 0.52917721092


@dataclass(frozen=True)
class PEDTable:
    """Potential-energy distribution in non-redundant GIC coordinates."""

    values: np.ndarray
    labels: tuple[str, ...]


@dataclass(frozen=True)
class InternalGFResult:
    frequencies_cm: np.ndarray
    force_constants: np.ndarray
    g_matrix: np.ndarray
    b_matrix: np.ndarray
    u_matrix: np.ndarray
    modes_internal: np.ndarray
    ped: PEDTable
    primitive_labels: tuple[str, ...]
    gic_labels: tuple[str, ...]
    gic_names: tuple[str, ...] = ()
    gic_irreps: tuple[str, ...] = ()
    point_group: str = "UNKNOWN"
    symmetrized_gics: bool = False
    scaling_factors: np.ndarray | None = None
    coordinate_source: str = "frozen-gic-definition"


def primitive_label(primitive: object) -> str:
    function = getattr(primitive, "function", "")
    atoms = tuple(int(atom) for atom in getattr(primitive, "atoms"))
    mode = int(getattr(primitive, "mode", 0))
    atoms_text = ",".join(str(atom) for atom in atoms)
    if function == "R":
        return f"R({atoms_text})"
    if function == "A":
        return f"A({atoms_text})"
    if function == "D":
        return f"D({atoms_text})"
    if function == "U":
        return f"U({atoms_text})"
    if function == "IMPD":
        center, n1, n2, n3 = atoms
        return f"D({n1},{center},{n3},{n2})"
    if function == "L":
        return f"L({atoms_text},0,{mode})"
    if function == "RPCK":
        return getattr(primitive, "gaussian_expression")()
    suffix = f":{mode}" if mode else ""
    return f"{function}{suffix}({atoms_text})"


def gic_labels_from_u(
    u_matrix: np.ndarray,
    primitive_labels: tuple[str, ...],
    *,
    threshold: float = 0.15,
) -> tuple[str, ...]:
    labels: list[str] = []
    for col in range(u_matrix.shape[1]):
        terms = []
        for row, coeff in enumerate(u_matrix[:, col]):
            if abs(coeff) < threshold:
                continue
            sign = "+" if coeff >= 0.0 else "-"
            terms.append(f"{sign}{abs(coeff):.3f}*{primitive_labels[row]}")
        labels.append(" ".join(terms) if terms else f"GIC{col + 1}")
    return tuple(labels)


def _mass_inverse(masses_amu: np.ndarray) -> np.ndarray:
    weights = np.repeat(1.0 / np.asarray(masses_amu, dtype=float), 3)
    return np.diag(weights)


def _internal_backtransform(bq: np.ndarray, masses_amu: np.ndarray, g_matrix: np.ndarray) -> np.ndarray:
    minv = _mass_inverse(masses_amu)
    return minv @ bq.T @ np.linalg.pinv(g_matrix, rcond=1.0e-10)


def _ped(force_constants: np.ndarray, modes_internal: np.ndarray, eigenvalues: np.ndarray) -> np.ndarray:
    ped = np.zeros((force_constants.shape[0], modes_internal.shape[1]), dtype=float)
    for mode in range(modes_internal.shape[1]):
        lam = eigenvalues[mode]
        if abs(lam) < 1.0e-14:
            continue
        vector = modes_internal[:, mode]
        raw = vector * (force_constants @ vector) / lam
        total = float(np.sum(np.abs(raw)))
        if total > 0.0:
            ped[:, mode] = 100.0 * np.abs(raw) / total
    return ped


def pulay_scale_internal_hessian(
    force_constants: np.ndarray,
    diagonal_factors: np.ndarray | None,
) -> np.ndarray:
    """Scale an internal-coordinate Hessian with Pulay-style factors."""
    f_mat = np.asarray(force_constants, dtype=float)
    if diagonal_factors is None:
        return np.array(f_mat, dtype=float, copy=True)
    factors = np.asarray(diagonal_factors, dtype=float)
    if factors.shape != (f_mat.shape[0],):
        raise ValueError(f"Scaling factors must have length {f_mat.shape[0]}")
    if np.any(factors < 0.0):
        raise ValueError("Pulay scaling factors must be non-negative")
    scale = np.sqrt(np.outer(factors, factors))
    return 0.5 * (f_mat * scale + (f_mat * scale).T)


def gf_from_cartesian_hessian_and_gic_b_matrix(
    cartesian_hessian: np.ndarray,
    b_matrix_internal: np.ndarray,
    masses_amu: np.ndarray,
    *,
    gic_labels: tuple[str, ...],
    primitive_labels: tuple[str, ...] = (),
    u_matrix: np.ndarray | None = None,
    gic_names: tuple[str, ...] = (),
    gic_irreps: tuple[str, ...] = (),
    point_group: str = "UNKNOWN",
    symmetrized_gics: bool = False,
    scaling_factors: np.ndarray | None = None,
    coordinate_source: str = "frozen-gic-definition",
) -> InternalGFResult:
    """Run Wilson GF from a Cartesian Hessian and a fixed non-redundant B matrix."""
    hessian = np.asarray(cartesian_hessian, dtype=float)
    bq = np.asarray(b_matrix_internal, dtype=float)
    masses = np.asarray(masses_amu, dtype=float)
    if hessian.shape != (3 * len(masses), 3 * len(masses)):
        raise ValueError("Cartesian Hessian has inconsistent dimensions")
    if bq.ndim != 2 or bq.shape[1] != hessian.shape[0]:
        raise ValueError("Internal B matrix has inconsistent dimensions")
    if len(gic_labels) != bq.shape[0]:
        raise ValueError("GIC label count does not match B matrix rows")

    minv = _mass_inverse(masses)
    g_matrix = bq @ minv @ bq.T
    backtransform = _internal_backtransform(bq, masses, g_matrix)
    force_constants = backtransform.T @ hessian @ backtransform
    force_constants = 0.5 * (force_constants + force_constants.T)
    force_constants = pulay_scale_internal_hessian(force_constants, scaling_factors)

    gf = solve_wilson_gf(force_constants, g_matrix, scale_to_cm=True)
    g_eval, g_vec = np.linalg.eigh(0.5 * (g_matrix + g_matrix.T))
    g_inv_half = (g_vec * (1.0 / np.sqrt(np.clip(g_eval, 1.0e-14, None)))) @ g_vec.T
    modes_internal = g_inv_half @ gf.normal_modes
    ped = _ped(force_constants, modes_internal, gf.eigenvalues)
    return InternalGFResult(
        frequencies_cm=gf.frequencies_cm,
        force_constants=force_constants,
        g_matrix=g_matrix,
        b_matrix=bq,
        u_matrix=np.array(
            u_matrix if u_matrix is not None else np.eye(bq.shape[0]),
            dtype=float,
            copy=True,
        ),
        modes_internal=modes_internal,
        ped=PEDTable(ped, gic_labels),
        primitive_labels=tuple(primitive_labels),
        gic_labels=tuple(gic_labels),
        gic_names=tuple(gic_names),
        gic_irreps=tuple(gic_irreps),
        point_group=point_group,
        symmetrized_gics=bool(symmetrized_gics),
        scaling_factors=None if scaling_factors is None else np.asarray(scaling_factors, dtype=float),
        coordinate_source=coordinate_source,
    )


def gf_from_hessian_input_and_gic_definition(
    input_data: HessianInput,
    definition: GICDefinition,
    *,
    coordinates_angstrom: np.ndarray | None = None,
    scaling_factors: np.ndarray | None = None,
) -> InternalGFResult:
    """Run GF/PED using a frozen GIC definition and canonical Hessian input."""
    input_data.validate()
    if len(definition.reference_coordinates_angstrom) != len(input_data.atomic_numbers):
        raise ValueError("GIC definition atom count does not match Hessian input")
    coords_for_b = (
        np.asarray(input_data.cartesian_coordinates_bohr, dtype=float)
        if coordinates_angstrom is None
        else np.asarray(coordinates_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    )
    b_matrix = build_gic_b_matrix(definition, coordinates_angstrom=coords_for_b)
    return gf_from_cartesian_hessian_and_gic_b_matrix(
        input_data.cartesian_hessian,
        np.asarray(b_matrix.rows, dtype=float),
        input_data.masses_amu,
        gic_labels=_gic_display_labels(definition),
        primitive_labels=tuple(primitive_label(primitive) for primitive in definition.primitives),
        gic_names=b_matrix.coordinate_names,
        gic_irreps=b_matrix.irreps,
        point_group=definition.point_group,
        symmetrized_gics=definition.symmetrize,
        scaling_factors=scaling_factors,
        coordinate_source=f"frozen-gic-definition:{definition.point_group}",
    )


def gf_from_hessian_input_and_xyzin(
    input_data: HessianInput,
    xyzin_path: Path,
    *,
    scaling_factors: np.ndarray | None = None,
) -> InternalGFResult:
    """Run GF/PED using the frozen #GIC section stored in an ORACLE xyzin file."""
    input_data.validate()
    definition = read_gic_definition_from_xyzin(Path(xyzin_path))
    b_matrix = build_gic_b_matrix(
        definition,
        coordinates_angstrom=np.asarray(input_data.cartesian_coordinates_bohr, dtype=float),
    )
    return gf_from_cartesian_hessian_and_gic_b_matrix(
        input_data.cartesian_hessian,
        np.asarray(b_matrix.rows, dtype=float),
        input_data.masses_amu,
        gic_labels=_gic_display_labels(definition),
        primitive_labels=tuple(primitive_label(primitive) for primitive in definition.primitives),
        gic_names=b_matrix.coordinate_names,
        gic_irreps=b_matrix.irreps,
        point_group=definition.point_group,
        symmetrized_gics=definition.symmetrize,
        scaling_factors=scaling_factors,
        coordinate_source=f"xyzin-frozen-gic:{Path(xyzin_path)}",
    )


def gf_from_cartesian_hessian_and_oracle_gics(
    cartesian_hessian: np.ndarray,
    coordinates_bohr: np.ndarray,
    atomic_numbers: np.ndarray,
    masses_amu: np.ndarray,
    *,
    symmetrize: bool = False,
) -> InternalGFResult:
    """Run Wilson GF from a Cartesian Hessian and freshly generated ORACLE GICs."""
    with tempfile.TemporaryDirectory(prefix="oracle-gf-") as tmp:
        tmpdir = Path(tmp)
        xyz = tmpdir / "geometry.xyz"
        xyzin = tmpdir / "geometry.xyzin"
        _write_xyz_from_hessian_geometry(xyz, atomic_numbers, coordinates_bohr)
        preprocess_to_enriched_xyz(xyz, xyzin)
        write_validation_section(xyzin)
        write_gicforge_build_sections(xyzin, symmetrize=symmetrize)
        result = gf_from_hessian_input_and_xyzin(
            HessianInput(
                atomic_numbers=np.asarray(atomic_numbers, dtype=int),
                cartesian_coordinates_bohr=np.asarray(coordinates_bohr, dtype=float),
                masses_amu=np.asarray(masses_amu, dtype=float),
                cartesian_hessian=np.asarray(cartesian_hessian, dtype=float),
                harmonic_frequencies_cm=np.array((), dtype=float),
                source="generated-oracle-gics",
            ),
            xyzin,
        )
    return InternalGFResult(
        frequencies_cm=result.frequencies_cm,
        force_constants=result.force_constants,
        g_matrix=result.g_matrix,
        b_matrix=result.b_matrix,
        u_matrix=result.u_matrix,
        modes_internal=result.modes_internal,
        ped=result.ped,
        primitive_labels=result.primitive_labels,
        gic_labels=result.gic_labels,
        gic_names=result.gic_names,
        gic_irreps=result.gic_irreps,
        point_group=result.point_group,
        symmetrized_gics=result.symmetrized_gics,
        scaling_factors=result.scaling_factors,
        coordinate_source="generated-oracle-gics",
    )


def gf_from_hessian_input_with_oracle_gics(input_data: HessianInput) -> InternalGFResult:
    """Run GF/PED from canonical Hessian input using freshly generated ORACLE GICs."""
    input_data.validate()
    return gf_from_cartesian_hessian_and_oracle_gics(
        input_data.cartesian_hessian,
        input_data.cartesian_coordinates_bohr,
        input_data.atomic_numbers,
        input_data.masses_amu,
    )


def gf_from_gaussian_fchk_with_oracle_gics(path: Path) -> InternalGFResult:
    """Gaussian adapter: read FCHK, then run the diagnostic ORACLE GF path."""
    from oracle_gaussian import hessian_input_from_gaussian_fchk

    return gf_from_hessian_input_with_oracle_gics(hessian_input_from_gaussian_fchk(path))


def _write_xyz_from_hessian_geometry(
    path: Path,
    atomic_numbers: np.ndarray,
    coordinates_bohr: np.ndarray,
) -> None:
    atoms = tuple(atomic_symbol(int(number)) for number in atomic_numbers)
    coords = np.asarray(coordinates_bohr, dtype=float) * BOHR_TO_ANGSTROM
    lines = [str(len(atoms)), "ORACLE GF generated geometry"]
    for atom, (x, y, z) in zip(atoms, coords):
        lines.append(f"{atom:2s} {x:16.10f} {y:16.10f} {z:16.10f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _gic_display_labels(definition: GICDefinition) -> tuple[str, ...]:
    labels: list[str] = []
    for gic in definition.gics:
        if gic.gaussian_expression and gic.gaussian_expression != "NONE":
            labels.append(f"{gic.identifier} {gic.gaussian_expression}")
        else:
            labels.append(gic.identifier)
    return tuple(labels)
