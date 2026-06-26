"""Molecular point-group detection and serialized symmetry operations."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import permutations, product
import re

import numpy as np

from .geometry import MolecularGeometry
from .topology.elements import atomic_number


@dataclass(frozen=True)
class SymmetryOperation:
    label: str
    rotation: tuple[tuple[float, float, float], ...]
    permutation: tuple[int, ...]
    max_deviation: float


@dataclass(frozen=True)
class MolecularSymmetry:
    point_group: str
    operations: tuple[SymmetryOperation, ...]
    atom_classes: tuple[tuple[int, ...], ...]
    max_deviation: float
    mean_deviation: float


def analyze_molecular_symmetry(
    geometry: MolecularGeometry,
    *,
    distance_tolerance: float,
    inertia_tolerance: float,
    max_rotation_order: int,
) -> MolecularSymmetry:
    symbols = list(geometry.atoms)
    weights = np.array([atomic_number(symbol) or 1 for symbol in symbols], dtype=float)
    oriented = orient_coords(geometry.coordinates_angstrom, weights=weights)
    elements, atom_classes, permutations = symmetry_elements_from_geometry(
        symbols,
        oriented,
        tol=distance_tolerance,
        max_n=max_rotation_order,
        tol_H=distance_tolerance,
        ignore_isotopes=True,
        auto_max_n=True,
        inertia_tol=inertia_tolerance,
    )
    if not elements:
        elements = [("E", np.eye(3), 0.0)]
        permutations = [tuple(range(len(symbols)))]
        atom_classes = tuple((idx,) for idx in range(len(symbols)))
    operations = tuple(
        SymmetryOperation(
            label=str(label),
            rotation=tuple(tuple(float(value) for value in row) for row in rotation),
            permutation=tuple(int(item) + 1 for item in permutation),
            max_deviation=float(max_deviation),
        )
        for (label, rotation, max_deviation), permutation in zip(elements, permutations)
    )
    return MolecularSymmetry(
        point_group=group_label(elements, linear=is_linear(oriented, tol=distance_tolerance)),
        operations=operations,
        atom_classes=tuple(tuple(int(atom) + 1 for atom in cls) for cls in atom_classes),
        max_deviation=float(max((op.max_deviation for op in operations), default=0.0)),
        mean_deviation=(
            float(np.mean([op.max_deviation for op in operations])) if operations else 0.0
        ),
    )


def symmetry_section_lines(symmetry: MolecularSymmetry, *, thresholds) -> list[str]:
    lines = [
        "SCHEMA oracle.xyz.symmetry.v1",
        f"POINT_GROUP {symmetry.point_group}",
        f"OPERATION_COUNT {len(symmetry.operations)}",
        f"MAX_OPERATION_DEVIATION_ANGSTROM {symmetry.max_deviation:.12g}",
        f"MEAN_OPERATION_DEVIATION_ANGSTROM {symmetry.mean_deviation:.12g}",
        f"THRESHOLD_DISTANCE_ANGSTROM {thresholds.distance_angstrom:.12g}",
        f"THRESHOLD_INERTIA_RELATIVE {thresholds.inertia_relative:.12g}",
        f"MAX_ROTATION_ORDER {thresholds.max_rotation_order}",
        "[OPERATIONS]",
    ]
    for idx, operation in enumerate(symmetry.operations, start=1):
        matrix = ",".join(
            f"{value:.12g}" for row in operation.rotation for value in row
        )
        permutation = ",".join(str(atom) for atom in operation.permutation)
        lines.append(
            f"{idx} LABEL={operation.label} "
            f"MAX_DEVIATION={operation.max_deviation:.12g} "
            f"PERMUTATION={permutation} MATRIX={matrix}"
        )
    lines.append("[ATOM_CLASSES]")
    if symmetry.atom_classes:
        for idx, atoms in enumerate(symmetry.atom_classes, start=1):
            lines.append(f"{idx} ATOMS=" + ",".join(str(atom) for atom in atoms))
    else:
        lines.append("NONE")
    return lines


def orient_coords(coords, weights=None):
    x = np.array(coords, dtype=float)
    w = np.ones(len(x), dtype=float) if weights is None else np.array(weights, dtype=float)
    center = np.sum(x * w[:, None], axis=0) / max(float(np.sum(w)), 1.0e-12)
    x = x - center
    inertia = np.zeros((3, 3), dtype=float)
    for idx, vec in enumerate(x):
        inertia += w[idx] * ((np.dot(vec, vec) * np.eye(3)) - np.outer(vec, vec))
    evals, evecs = np.linalg.eigh(inertia)
    frame = evecs[:, np.argsort(evals)]
    if np.linalg.det(frame) < 0.0:
        frame[:, -1] *= -1.0
    return x @ frame


def is_linear(coords, tol=1.0e-3):
    x = np.array(coords, dtype=float)
    inertia = np.zeros((3, 3), dtype=float)
    for vec in x:
        inertia += np.dot(vec, vec) * np.eye(3) - np.outer(vec, vec)
    return bool(np.linalg.eigvalsh(inertia)[0] < tol)


def symmetry_elements_from_geometry(
    symbols,
    coords_oriented,
    tol=1.0e-3,
    max_n=6,
    tol_H=None,
    ignore_isotopes=False,
    auto_max_n=False,
    inertia_tol=1.0e-3,
):
    coords = np.asarray(coords_oriented, dtype=float)
    radii = np.linalg.norm(coords, axis=1)
    max_radius = float(np.max(radii)) if len(radii) else 1.0
    if max_radius <= 0.0:
        max_radius = 1.0
    scaled = coords / max_radius
    sym_use = [symbol[0] for symbol in symbols] if ignore_isotopes else list(symbols)
    if auto_max_n:
        inertia = np.zeros((3, 3), dtype=float)
        for vec in coords:
            inertia += np.dot(vec, vec) * np.eye(3) - np.outer(vec, vec)
        evals = np.linalg.eigvalsh(inertia)
        max_inertia = float(np.max(evals)) if len(evals) else 0.0
        if max_inertia > 0.0:
            d01 = abs(evals[0] - evals[1]) / max_inertia
            d12 = abs(evals[1] - evals[2]) / max_inertia
            if d01 > inertia_tol and d12 > inertia_tol:
                max_n = min(max_n, 2)
    elements = []
    permutations = []
    seen: set[tuple[tuple[int, ...], tuple[float, ...]]] = set()
    for label, rotation in candidate_ops(max_n=max_n):
        mapped, max_dev = _match_with_map(
            sym_use,
            scaled,
            scaled @ rotation.T,
            tol,
            tol_H=tol_H,
        )
        if mapped is not None:
            unique_key = (
                tuple(int(item) for item in mapped),
                tuple(round(float(value), 10) for value in rotation.reshape(-1)),
            )
            if unique_key in seen:
                continue
            seen.add(unique_key)
            elements.append((label, rotation, float(max_dev)))
            permutations.append(tuple(mapped))
    return elements, _atom_classes(len(symbols), permutations), permutations


def group_label(elements, linear=False):
    labels = [item[0] for item in elements]
    polyhedral = _polyhedral_group_label(elements)
    if polyhedral:
        return polyhedral
    dnd_n = _dnd_group_order(labels)
    if dnd_n:
        return f"D{dnd_n}d"
    nmax, axis = _highest_cn_axis(labels)
    has_i = "i" in labels
    has_sigma = any(label.startswith("sigma") for label in labels)
    has_c2 = any(label.startswith("C2") for label in labels)
    if linear:
        return "Dinfh" if has_i else "Cinfv"
    if nmax >= 2:
        sigma_h = {"x": "sigma_yz", "y": "sigma_xz", "z": "sigma_xy"}.get(axis or "z")
        has_sigma_h = sigma_h in labels
        has_sigma_v = has_sigma and not has_sigma_h
        if has_sigma_h and has_c2:
            return f"D{nmax}h"
        if has_sigma_h:
            return f"C{nmax}h"
        if has_sigma_v:
            return f"C{nmax}v"
        if has_c2:
            return f"D{nmax}"
        return f"C{nmax}"
    if has_i:
        return "Ci"
    if has_sigma:
        return "Cs"
    return "C1"


def _match_with_map(symbols, coords1, coords2, tol, tol_H=None):
    used = np.zeros(len(coords2), dtype=bool)
    mapping = [-1] * len(coords1)
    by_symbol: dict[str, list[int]] = {}
    for idx, symbol in enumerate(symbols):
        by_symbol.setdefault(symbol, []).append(idx)
    radii2 = np.linalg.norm(coords2, axis=1)
    max_dev = 0.0
    for idx in sorted(range(len(coords1)), key=lambda item: (len(by_symbol[symbols[item]]), item)):
        eff_tol = tol_H if tol_H is not None and symbols[idx] == "H" else tol
        radius = float(np.linalg.norm(coords1[idx]))
        candidates = [
            cand
            for cand in by_symbol.get(symbols[idx], ())
            if not used[cand] and abs(radii2[cand] - radius) <= eff_tol
        ]
        candidates.sort(key=lambda cand: abs(radii2[cand] - radius))
        for cand in candidates:
            deviation = float(np.linalg.norm(coords1[idx] - coords2[cand]))
            if deviation < eff_tol:
                mapping[idx] = cand
                used[cand] = True
                max_dev = max(max_dev, deviation)
                break
        if mapping[idx] < 0:
            return None, None
    return tuple(mapping), max_dev


def _atom_classes(natoms: int, permutations) -> tuple[tuple[int, ...], ...]:
    parent = list(range(natoms))

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left, right):
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for permutation in permutations:
        for left, right in enumerate(permutation):
            union(left, int(right))
    classes: dict[int, list[int]] = {}
    for atom in range(natoms):
        classes.setdefault(find(atom), []).append(atom)
    return tuple(tuple(values) for values in classes.values())


def _highest_cn_axis(labels):
    best = 1
    axis = None
    for label in labels:
        match = re.match(r"C(\d+)([xyz])", str(label))
        if match and int(match.group(1)) > best:
            best = int(match.group(1))
            axis = match.group(2)
    return best, axis


def _polyhedral_group_label(elements) -> str | None:
    matrices = [np.asarray(item[1], dtype=float) for item in elements]
    td_classes = Counter(
        operation_class
        for operation_class in (_td_candidate_class(matrix) for matrix in matrices)
        if operation_class is not None
    )
    if td_classes == Counter(
        {
            "E": 1,
            "C3": 8,
            "C2": 3,
            "S4": 6,
            "sigma_d": 6,
        }
    ):
        return "Td"
    oh_classes = Counter(
        operation_class
        for operation_class in (_polyhedral_matrix_class(matrix) for matrix in matrices)
        if operation_class is not None
    )
    if oh_classes == Counter(
        {
            "E": 1,
            "C3": 8,
            "C2_axis": 3,
            "C4": 6,
            "C2_edge": 6,
            "i*E": 1,
            "i*C3": 8,
            "i*C2_axis": 3,
            "i*C4": 6,
            "i*C2_edge": 6,
        }
    ):
        return "Oh"
    ih_classes = Counter(
        operation_class
        for operation_class in (_icosahedral_candidate_class(matrix) for matrix in matrices)
        if operation_class is not None
    )
    if ih_classes == Counter(
        {
            "E": 1,
            "C2": 15,
            "C3": 20,
            "C5": 12,
            "C5_2": 12,
            "i*E": 1,
            "i*C2": 15,
            "i*C3": 20,
            "i*C5": 12,
            "i*C5_2": 12,
        }
    ):
        return "Ih"
    if ih_classes == Counter({"E": 1, "C2": 15, "C3": 20, "C5": 12, "C5_2": 12}):
        return "I"
    return None


def _polyhedral_matrix_class(matrix: np.ndarray) -> str | None:
    if np.allclose(matrix, np.eye(3), atol=1.0e-8):
        return "E"
    det = float(np.linalg.det(matrix))
    trace = float(np.trace(matrix))
    if det > 0.0:
        if abs(trace) <= 1.0e-8:
            return "C3"
        if abs(trace - 1.0) <= 1.0e-8:
            return "C4"
        if abs(trace + 1.0) <= 1.0e-8:
            return "C2_axis" if _is_coordinate_axis_c2(matrix) else "C2_edge"
    if det < 0.0:
        if np.allclose(matrix, -np.eye(3), atol=1.0e-8):
            return "i*E"
        proper = -matrix
        proper_class = _polyhedral_matrix_class(proper)
        if proper_class is not None:
            return f"i*{proper_class}"
        if abs(trace - 1.0) <= 1.0e-8:
            return "sigma_d"
        if abs(trace + 1.0) <= 1.0e-8:
            return "S4"
    return None


def _is_coordinate_axis_c2(matrix: np.ndarray) -> bool:
    rounded = np.rint(matrix)
    if not np.allclose(matrix, rounded, atol=1.0e-8):
        return False
    return bool(np.count_nonzero(np.abs(np.diag(rounded)) > 0.5) == 3)


def _dnd_group_order(labels: list[str]) -> int | None:
    orders = []
    for label in labels:
        for pattern in (r"Dnd_C(\d+)z\^", r"Dnd_C2_xy_(\d+)_"):
            match = re.search(pattern, str(label))
            if match:
                orders.append(int(match.group(1)) // 2)
    return max(orders) if orders else None


@lru_cache(maxsize=16)
def candidate_ops(max_n=6):
    ops = [("E", np.eye(3)), ("i", -np.eye(3))]
    for axis, name in [(0, "sigma_yz"), (1, "sigma_xz"), (2, "sigma_xy")]:
        matrix = np.eye(3)
        matrix[axis, axis] = -1.0
        ops.append((name, matrix))
    for n in range(2, max_n + 1):
        for power in range(1, n):
            theta = 2.0 * np.pi * power / n
            ops.append((f"C{n}z^{power}", _rotation_matrix((0, 0, 1), theta)))
            ops.append((f"C{n}x^{power}", _rotation_matrix((1, 0, 0), theta)))
            ops.append((f"C{n}y^{power}", _rotation_matrix((0, 1, 0), theta)))
            sigma_h = np.diag((1.0, 1.0, -1.0))
            ops.append(
                (
                    f"sigma_h*C{n}z^{power}",
                    sigma_h @ _rotation_matrix((0, 0, 1), theta),
                )
            )
    for n in range(2, max_n + 1):
        for k in range(n):
            theta = np.pi * k / n
            ops.append(
                (
                    f"C2_xy_{n}_{k}",
                    _rotation_matrix((np.cos(theta), np.sin(theta), 0), np.pi),
                )
            )
    for n in range(3, max_n + 1):
        for k in range(n):
            theta = np.pi * k / n
            normal = np.array((np.cos(theta), np.sin(theta), 0.0), dtype=float)
            normal /= np.linalg.norm(normal)
            ops.append((f"sigma_v_{n}_{k}", np.eye(3) - 2.0 * np.outer(normal, normal)))
    ops.extend(_dnd_candidate_ops(max_n))
    ops.extend(_cubic_candidate_ops())
    ops.extend(_icosahedral_candidate_ops())
    return ops


def _dnd_candidate_ops(max_n: int) -> list[tuple[str, np.ndarray]]:
    ops: list[tuple[str, np.ndarray]] = []
    sd = _diagonal_reflection_matrix()
    for n in range(2, max_n + 1):
        order = 2 * n
        for power in range(n):
            theta = 2.0 * np.pi * power / n
            ops.append(
                (
                    f"Dnd_C2_xy_{order}_{(2 * power + 1) % order}",
                    sd @ _rotation_matrix((0, 0, 1), theta),
                )
            )
        for k in range(n):
            theta = np.pi * k / n
            c2 = _rotation_matrix((np.cos(theta), np.sin(theta), 0), np.pi)
            ops.append((f"Dnd_C{order}z^{(2 * k + 1) % order}", sd @ c2))
    return ops


def _diagonal_reflection_matrix() -> np.ndarray:
    return np.array(
        ((0.0, 1.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        dtype=float,
    )


def _cubic_candidate_ops() -> list[tuple[str, np.ndarray]]:
    ops: list[tuple[str, np.ndarray]] = []
    counters: Counter[str] = Counter()
    seen: set[tuple[float, ...]] = set()

    for matrix in _signed_permutation_matrices():
        if _sign_product(matrix) != 1:
            continue
        label = _unique_cubic_label(_td_candidate_class(matrix), "td", counters)
        if label is not None:
            key = tuple(float(value) for value in matrix.reshape(-1))
            seen.add(key)
            ops.append((label, matrix))

    for matrix in _signed_permutation_matrices():
        key = tuple(float(value) for value in matrix.reshape(-1))
        if key in seen:
            continue
        label = _unique_cubic_label(_oh_candidate_class(matrix), "oh", counters)
        if label is not None:
            ops.append((label, matrix))
    return ops


def _signed_permutation_matrices() -> tuple[np.ndarray, ...]:
    matrices: list[np.ndarray] = []
    for permutation in permutations(range(3)):
        for signs in product((-1.0, 1.0), repeat=3):
            matrix = np.zeros((3, 3), dtype=float)
            for row, column in enumerate(permutation):
                matrix[row, column] = signs[row]
            matrices.append(matrix)
    return tuple(matrices)


def _sign_product(matrix: np.ndarray) -> int:
    product_value = 1.0
    for row in range(3):
        nonzero = np.flatnonzero(np.abs(matrix[row]) > 0.5)
        if len(nonzero) != 1:
            return 0
        product_value *= float(matrix[row, nonzero[0]])
    return 1 if product_value > 0.0 else -1


def _td_candidate_class(matrix: np.ndarray) -> str | None:
    if np.allclose(matrix, np.eye(3), atol=1.0e-8):
        return "E"
    det = float(np.linalg.det(matrix))
    trace = float(np.trace(matrix))
    if det > 0.0:
        if abs(trace) <= 1.0e-8:
            return "C3"
        if abs(trace + 1.0) <= 1.0e-8:
            return "C2"
    if det < 0.0:
        if abs(trace + 1.0) <= 1.0e-8:
            return "S4"
        if abs(trace - 1.0) <= 1.0e-8:
            return "sigma_d"
    return None


def _oh_candidate_class(matrix: np.ndarray) -> str | None:
    operation_class = _polyhedral_matrix_class(matrix)
    if operation_class is None:
        return None
    return operation_class.replace("*", "_")


def _unique_cubic_label(
    operation_class: str | None,
    prefix: str,
    counters: Counter[str],
) -> str | None:
    if operation_class is None:
        return None
    if operation_class == "E":
        return f"{prefix}_E"
    if operation_class == "i_E":
        return "i"
    counters[f"{prefix}_{operation_class}"] += 1
    return f"{prefix}_{operation_class}_{counters[f'{prefix}_{operation_class}']}"


def _icosahedral_candidate_ops() -> list[tuple[str, np.ndarray]]:
    ops: list[tuple[str, np.ndarray]] = []
    counters: Counter[str] = Counter()
    for matrix in _icosahedral_operation_matrices():
        operation_class = _icosahedral_candidate_class(matrix)
        if operation_class is None:
            continue
        label = _unique_icosahedral_label(operation_class, counters)
        ops.append((label, matrix))
    return ops


@lru_cache(maxsize=1)
def _icosahedral_operation_matrices() -> tuple[np.ndarray, ...]:
    vertices = _icosahedral_vertices()
    source_indices = _independent_vertex_triple(vertices)
    source = vertices[list(source_indices)]
    source_inverse = np.linalg.inv(source)
    matrices: list[np.ndarray] = []
    seen: set[tuple[float, ...]] = set()
    for target_indices in permutations(range(len(vertices)), 3):
        target = vertices[list(target_indices)]
        if abs(float(np.linalg.det(target))) <= 1.0e-8:
            continue
        matrix = target.T @ source_inverse.T
        if not np.allclose(matrix.T @ matrix, np.eye(3), atol=1.0e-8):
            continue
        if not _maps_vertices_to_self(matrix, vertices):
            continue
        key = tuple(round(float(value), 10) for value in matrix.reshape(-1))
        if key in seen:
            continue
        seen.add(key)
        matrices.append(matrix)
    return tuple(matrices)


def _icosahedral_vertices() -> np.ndarray:
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    vertices = []
    for y in (-1.0, 1.0):
        for z in (-phi, phi):
            vertices.append((0.0, y, z))
    for x in (-1.0, 1.0):
        for y in (-phi, phi):
            vertices.append((x, y, 0.0))
    for x in (-phi, phi):
        for z in (-1.0, 1.0):
            vertices.append((x, 0.0, z))
    return np.array(vertices, dtype=float)


def _independent_vertex_triple(vertices: np.ndarray) -> tuple[int, int, int]:
    for candidate in permutations(range(len(vertices)), 3):
        if abs(float(np.linalg.det(vertices[list(candidate)]))) > 1.0e-8:
            return tuple(int(item) for item in candidate)
    raise ValueError("no independent icosahedral vertex triple")


def _maps_vertices_to_self(matrix: np.ndarray, vertices: np.ndarray) -> bool:
    transformed = vertices @ matrix.T
    for vertex in transformed:
        if not np.any(np.all(np.isclose(vertices, vertex, atol=1.0e-8), axis=1)):
            return False
    return True


def _icosahedral_candidate_class(matrix: np.ndarray) -> str | None:
    if np.allclose(matrix, np.eye(3), atol=1.0e-8):
        return "E"
    det = float(np.linalg.det(matrix))
    proper = matrix if det > 0.0 else -matrix
    proper_class = _icosahedral_proper_class(proper)
    if proper_class is None:
        return None
    return f"i*{proper_class}" if det < 0.0 else proper_class


def _icosahedral_proper_class(matrix: np.ndarray) -> str | None:
    if np.allclose(matrix, np.eye(3), atol=1.0e-8):
        return "E"
    trace = float(np.trace(matrix))
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    phi_bar = (1.0 - np.sqrt(5.0)) / 2.0
    if abs(trace) <= 1.0e-8:
        return "C3"
    if abs(trace + 1.0) <= 1.0e-8:
        return "C2"
    if abs(trace - phi) <= 1.0e-8:
        return "C5"
    if abs(trace - phi_bar) <= 1.0e-8:
        return "C5_2"
    return None


def _unique_icosahedral_label(operation_class: str, counters: Counter[str]) -> str:
    if operation_class == "E":
        return "ih_E"
    if operation_class == "i*E":
        return "i"
    key = operation_class.replace("*", "_")
    counters[key] += 1
    return f"ih_{key}_{counters[key]}"


def _rotation_matrix(axis, theta):
    axis = np.array(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(theta)
    s = np.sin(theta)
    one_c = 1.0 - c
    return np.array(
        [
            [c + x * x * one_c, x * y * one_c - z * s, x * z * one_c + y * s],
            [y * x * one_c + z * s, c + y * y * one_c, y * z * one_c - x * s],
            [z * x * one_c - y * s, z * y * one_c + x * s, c + z * z * one_c],
        ],
        dtype=float,
    )
