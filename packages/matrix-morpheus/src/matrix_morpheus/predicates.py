from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from matrix_chem.topology.elements import atomic_number
from matrix_chem.topology.pipeline import build_topology_objects
from matrix_neo.survibfit.primitives import Primitive, eval_primitives

from .contracts import QMParameterPredicate


DEFAULT_INITIAL_GEOMETRY_PREDICATE_SCOPE = (
    "heavy_bonds",
    "oh_bonds",
    "heavy_angles",
    "coh_angles",
    "ring_torsions",
)


def initial_geometry_predicates(
    atoms: tuple[str, ...] | list[str],
    coords_angstrom: np.ndarray,
    *,
    distance_sigma_angstrom: float,
    angle_sigma_degree: float,
    dihedral_sigma_degree: float,
    scope: Iterable[str] = DEFAULT_INITIAL_GEOMETRY_PREDICATE_SCOPE,
    source: str = "initial_geometry",
) -> tuple[QMParameterPredicate, ...]:
    """Build primitive predicates centered on the initial Cartesian geometry."""

    if distance_sigma_angstrom <= 0.0:
        raise ValueError("distance_sigma_angstrom must be positive")
    if angle_sigma_degree <= 0.0:
        raise ValueError("angle_sigma_degree must be positive")
    if dihedral_sigma_degree <= 0.0:
        raise ValueError("dihedral_sigma_degree must be positive")
    symbols = tuple(str(atom) for atom in atoms)
    coords = np.asarray(coords_angstrom, dtype=float)
    z_numbers = np.array([atomic_number(symbol) or 0 for symbol in symbols], dtype=int)
    _continuous, graph, ringset, _synthons, _aromaticity = build_topology_objects(coords, z_numbers)
    scopes = {str(item).strip().lower() for item in scope if str(item).strip()}
    if "all" in scopes:
        scopes = set(DEFAULT_INITIAL_GEOMETRY_PREDICATE_SCOPE)
    predicates: list[QMParameterPredicate] = []

    bonds = tuple(sorted(tuple(sorted((int(i), int(j)))) for i, j in graph.bonds))
    if {"heavy_bonds", "oh_bonds"} & scopes:
        for i, j in bonds:
            zi, zj = int(z_numbers[i]), int(z_numbers[j])
            include = ("heavy_bonds" in scopes and zi > 1 and zj > 1) or (
                "oh_bonds" in scopes and {zi, zj} == {1, 8}
            )
            if not include:
                continue
            value = float(eval_primitives([Primitive("bond", (i, j))], coords)[0])
            predicates.append(
                QMParameterPredicate(
                    f"R({i + 1},{j + 1})",
                    value,
                    distance_sigma_angstrom,
                    f"{source}_distance",
                )
            )

    adjacency = tuple(
        tuple(sorted(int(item) for item in graph.adjacency[index]))
        for index in range(len(symbols))
    )
    if {"heavy_angles", "coh_angles"} & scopes:
        for center, neighbors in enumerate(adjacency):
            for pos, left in enumerate(neighbors):
                for right in neighbors[pos + 1 :]:
                    z_triplet = (int(z_numbers[left]), int(z_numbers[center]), int(z_numbers[right]))
                    include = (
                        "heavy_angles" in scopes and all(z > 1 for z in z_triplet)
                    ) or (
                        "coh_angles" in scopes
                        and z_triplet[1] == 8
                        and {z_triplet[0], z_triplet[2]} == {1, 6}
                    )
                    if not include:
                        continue
                    primitive = Primitive("angle", (left, center, right))
                    value = float(np.rad2deg(eval_primitives([primitive], coords)[0]))
                    predicates.append(
                        QMParameterPredicate(
                            f"A({left + 1},{center + 1},{right + 1})",
                            value,
                            angle_sigma_degree,
                            f"{source}_angle",
                        )
                    )

    if "ring_torsions" in scopes:
        for ring in getattr(ringset, "rings", ()):
            ring_atoms = tuple(int(atom) for atom in getattr(ring, "atoms", ()))
            size = len(ring_atoms)
            if size < 4:
                continue
            for idx in range(size):
                torsion = tuple(ring_atoms[(idx + offset) % size] for offset in range(4))
                value = float(np.rad2deg(eval_primitives([Primitive("dihedral", torsion)], coords)[0]))
                predicates.append(
                    QMParameterPredicate(
                        "D(" + ",".join(str(atom + 1) for atom in torsion) + ")",
                        value,
                        dihedral_sigma_degree,
                        f"{source}_ring_torsion",
                    )
                )

    return _unique_predicates(tuple(predicates))


def _unique_predicates(predicates: tuple[QMParameterPredicate, ...]) -> tuple[QMParameterPredicate, ...]:
    result: list[QMParameterPredicate] = []
    seen: set[tuple[str, str]] = set()
    for predicate in predicates:
        key = (predicate.label_pattern, predicate.source)
        if key in seen:
            continue
        seen.add(key)
        result.append(predicate)
    return tuple(result)
