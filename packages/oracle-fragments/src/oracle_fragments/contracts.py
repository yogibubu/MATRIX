from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from oracle_chem import read_enriched_xyz
from oracle_core import read_sectioned_lines, replace_section, section_content


ORACLE_XYZ_FRAGMENTS_SCHEMA = "oracle.xyz.fragments.v1"
ORACLE_XYZ_FRAGMENT_LIBRARY_SCHEMA = "oracle.xyz.fragment_library.v1"
ORACLE_XYZ_ASSEMBLY_SCHEMA = "oracle.xyz.assembly.v1"
ORACLE_XYZ_INTERACTION_CENTERS_SCHEMA = "oracle.xyz.interaction_centers.v1"

REQUIRED_TOPOLOGY_SCHEMA = "oracle.xyz.topology.v1"
REQUIRED_SYNTHONS_SCHEMA = "oracle.xyz.synthons.v1"
RANK_TOLERANCE = 1.0e-8


class FragmentContractError(ValueError):
    """Raised when a file cannot enter the ORACLE fragment workflow."""


@dataclass(frozen=True)
class FragmentRecord:
    identifier: str
    label: str
    atoms: tuple[int, ...]
    center: tuple[float, float, float]
    frame: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True)
class FragmentDefinition:
    strategy: str
    reference_fragment: str
    fragments: tuple[FragmentRecord, ...]


@dataclass(frozen=True)
class InteractionCenterRecord:
    identifier: str
    kind: str
    label: str
    atoms: tuple[int, ...]
    center: tuple[float, float, float]
    source: str


@dataclass(frozen=True)
class AtomCenterInteractionRecord:
    identifier: str
    kind: str
    atom: int
    center_id: str
    score: float
    source: str


@dataclass(frozen=True)
class InteractionCenterDefinition:
    strategy: str
    centers: tuple[InteractionCenterRecord, ...]
    interactions: tuple[AtomCenterInteractionRecord, ...]


def validate_fragment_prerequisites(path: Path) -> None:
    """Require saved topology and synthons before any fragment workflow starts."""
    lines = read_sectioned_lines(Path(path))
    _require_schema(lines, "TOPOLOGY", REQUIRED_TOPOLOGY_SCHEMA)
    _require_schema(lines, "SYNTHONS", REQUIRED_SYNTHONS_SCHEMA)


def fragment_plan_section_lines(
    *,
    status: str = "PLANNED",
    strategy: str = "TOPOLOGY_SYNTHON",
) -> list[str]:
    """Return the initial #FRAGMENTS section without computing fragments yet."""
    return [
        f"SCHEMA {ORACLE_XYZ_FRAGMENTS_SCHEMA}",
        f"STATUS {status.strip().upper()}",
        "DEPENDENCIES TOPOLOGY=oracle.xyz.topology.v1 SYNTHONS=oracle.xyz.synthons.v1",
        "INDEXING ATOMS=ONE_BASED",
        f"STRATEGY {strategy.strip().upper()}",
        "[FRAGMENTS]",
        "PENDING ROBUST_TOPOLOGY_CONTRACT",
    ]


def write_fragment_plan_section(path: Path) -> None:
    """Mark an enriched XYZ as ready for future topology-backed fragmentation."""
    target = Path(path)
    validate_fragment_prerequisites(target)
    replace_section(target, "FRAGMENTS", fragment_plan_section_lines())


def build_fragment_definition_from_xyzin(path: Path) -> FragmentDefinition:
    """Build concrete fragments from the saved topology connected components."""
    target = Path(path)
    validate_fragment_prerequisites(target)
    lines = read_sectioned_lines(target)
    geometry = read_enriched_xyz(target)
    bonds = _topology_bonds(lines, natoms=geometry.natoms)
    components = _connected_components(bonds, natoms=geometry.natoms)
    coords = np.asarray(geometry.coordinates_angstrom, dtype=float)
    fragments = tuple(
        FragmentRecord(
            identifier=f"F{idx:03d}",
            label=f"component_{idx}",
            atoms=tuple(component),
            center=_center(coords, component),
            frame=_frame(coords, component),
        )
        for idx, component in enumerate(components, start=1)
    )
    reference = max(fragments, key=lambda item: (len(item.atoms), -int(item.identifier[1:])))
    return FragmentDefinition(
        strategy="CONNECTED_COMPONENTS",
        reference_fragment=reference.identifier,
        fragments=fragments,
    )


def fragment_build_section_lines(definition: FragmentDefinition) -> list[str]:
    lines = [
        f"SCHEMA {ORACLE_XYZ_FRAGMENTS_SCHEMA}",
        "STATUS BUILT",
        "DEPENDENCIES TOPOLOGY=oracle.xyz.topology.v1 SYNTHONS=oracle.xyz.synthons.v1",
        "INDEXING ATOMS=ONE_BASED",
        f"STRATEGY {definition.strategy}",
        f"FRAGMENT_COUNT {len(definition.fragments)}",
        f"REFERENCE_FRAGMENT {definition.reference_fragment}",
        "[FRAGMENTS]",
    ]
    for fragment in definition.fragments:
        atoms = ",".join(str(atom) for atom in fragment.atoms)
        lines.append(
            f"{fragment.identifier} LABEL={fragment.label} SIZE={len(fragment.atoms)} ATOMS={atoms}"
        )
    lines.append("[CENTERS]")
    for fragment in definition.fragments:
        x, y, z = fragment.center
        lines.append(f"{fragment.identifier} X={x:.12g} Y={y:.12g} Z={z:.12g}")
    lines.append("[FRAMES]")
    for fragment in definition.fragments:
        axes = []
        for label, axis in zip(("X", "Y", "Z"), fragment.frame):
            axes.append(f"{label}={axis[0]:.12g},{axis[1]:.12g},{axis[2]:.12g}")
        lines.append(f"{fragment.identifier} {' '.join(axes)}")
    return lines


def write_fragment_build_section(path: Path) -> FragmentDefinition:
    """Materialize the #FRAGMENTS section from saved topology."""
    target = Path(path)
    definition = build_fragment_definition_from_xyzin(target)
    replace_section(target, "FRAGMENTS", fragment_build_section_lines(definition))
    return definition


def build_interaction_center_definition_from_xyzin(path: Path) -> InteractionCenterDefinition:
    """Build topology-backed bond/ring centers and atom-center candidates."""
    target = Path(path)
    validate_fragment_prerequisites(target)
    lines = read_sectioned_lines(target)
    geometry = read_enriched_xyz(target)
    coords = np.asarray(geometry.coordinates_angstrom, dtype=float)
    bonds = _topology_bonds(lines, natoms=geometry.natoms)
    rings = _topology_rings(lines, natoms=geometry.natoms)
    centers: list[InteractionCenterRecord] = []

    for left, right in bonds:
        centers.append(
            InteractionCenterRecord(
                identifier=f"C{len(centers) + 1:03d}",
                kind="BOND_CENTER",
                label=f"bond_{left}_{right}",
                atoms=(left, right),
                center=_center(coords, (left, right)),
                source="TOPOLOGY_BOND",
            )
        )
    for ring_index, atoms in rings:
        centers.append(
            InteractionCenterRecord(
                identifier=f"C{len(centers) + 1:03d}",
                kind="RING_CENTER",
                label=f"ring_{ring_index}",
                atoms=atoms,
                center=_center(coords, atoms),
                source="TOPOLOGY_RING",
            )
        )

    interactions = _atom_center_interactions(
        tuple(centers),
        bonds=bonds,
        coords=coords,
        natoms=geometry.natoms,
    )
    return InteractionCenterDefinition(
        strategy="TOPOLOGY_EQUIDISTANT_CENTER_CANDIDATES",
        centers=tuple(centers),
        interactions=interactions,
    )


def interaction_center_section_lines(definition: InteractionCenterDefinition) -> list[str]:
    lines = [
        f"SCHEMA {ORACLE_XYZ_INTERACTION_CENTERS_SCHEMA}",
        "STATUS BUILT",
        "DEPENDENCIES TOPOLOGY=oracle.xyz.topology.v1 SYNTHONS=oracle.xyz.synthons.v1",
        "INDEXING ATOMS=ONE_BASED",
        f"STRATEGY {definition.strategy}",
        f"CENTER_COUNT {len(definition.centers)}",
        f"INTERACTION_COUNT {len(definition.interactions)}",
        "[CENTERS]",
    ]
    if definition.centers:
        for center in definition.centers:
            atoms = ",".join(str(atom) for atom in center.atoms)
            x, y, z = center.center
            lines.append(
                f"{center.identifier} KIND={center.kind} LABEL={center.label} "
                f"ATOMS={atoms} X={x:.12g} Y={y:.12g} Z={z:.12g} SOURCE={center.source}"
            )
    else:
        lines.append("NONE")
    lines.append("[INTERACTIONS]")
    if definition.interactions:
        for interaction in definition.interactions:
            lines.append(
                f"{interaction.identifier} KIND={interaction.kind} ATOM={interaction.atom} "
                f"CENTER={interaction.center_id} SCORE={interaction.score:.8g} "
                f"SOURCE={interaction.source}"
            )
    else:
        lines.append("NONE")
    return lines


def write_interaction_center_section(path: Path) -> InteractionCenterDefinition:
    """Materialize virtual bond/ring centers and atom-center interaction candidates."""
    target = Path(path)
    definition = build_interaction_center_definition_from_xyzin(target)
    replace_section(target, "INTERACTION_CENTERS", interaction_center_section_lines(definition))
    return definition


def read_fragment_records(path: Path) -> tuple[FragmentRecord, ...]:
    """Read built fragment records from an enriched XYZ file."""
    lines = read_sectioned_lines(Path(path))
    section = section_content(lines, "FRAGMENTS")
    if not section:
        return ()
    if section[0].strip() != f"SCHEMA {ORACLE_XYZ_FRAGMENTS_SCHEMA}":
        raise FragmentContractError("invalid #FRAGMENTS schema")
    status = _section_value(section, "STATUS")
    if status != "BUILT":
        return ()
    fragment_rows = _subsection(section, "FRAGMENTS")
    centers = _center_rows(_subsection(section, "CENTERS"))
    frames = _frame_rows(_subsection(section, "FRAMES"))
    records: list[FragmentRecord] = []
    for row in fragment_rows:
        parts = row.split()
        if not parts:
            continue
        identifier = parts[0]
        fields = _key_values(parts[1:])
        atoms_text = fields.get("ATOMS", "")
        if not atoms_text:
            raise FragmentContractError(f"fragment {identifier} has no ATOMS field")
        try:
            atoms = tuple(int(item) for item in atoms_text.split(",") if item)
        except ValueError as exc:
            raise FragmentContractError(f"invalid ATOMS field for {identifier}") from exc
        records.append(
            FragmentRecord(
                identifier=identifier,
                label=fields.get("LABEL", identifier),
                atoms=atoms,
                center=centers.get(identifier, (0.0, 0.0, 0.0)),
                frame=frames.get(identifier, _identity_frame()),
            )
        )
    return tuple(records)


def read_interaction_center_definition(path: Path) -> InteractionCenterDefinition:
    """Read a built #INTERACTION_CENTERS section, if present."""
    lines = read_sectioned_lines(Path(path))
    section = section_content(lines, "INTERACTION_CENTERS")
    if not section:
        return InteractionCenterDefinition(strategy="NONE", centers=(), interactions=())
    if section[0].strip() != f"SCHEMA {ORACLE_XYZ_INTERACTION_CENTERS_SCHEMA}":
        raise FragmentContractError("invalid #INTERACTION_CENTERS schema")
    status = _section_value(section, "STATUS")
    if status != "BUILT":
        return InteractionCenterDefinition(strategy="NONE", centers=(), interactions=())
    centers: list[InteractionCenterRecord] = []
    for row in _subsection(section, "CENTERS"):
        if row.strip().upper() == "NONE":
            continue
        parts = row.split()
        fields = _key_values(parts[1:])
        atoms = _parse_int_list(fields.get("ATOMS", ""))
        try:
            centers.append(
                InteractionCenterRecord(
                    identifier=parts[0],
                    kind=fields.get("KIND", "UNKNOWN"),
                    label=fields.get("LABEL", parts[0]),
                    atoms=atoms,
                    center=(
                        float(fields.get("X", "0.0")),
                        float(fields.get("Y", "0.0")),
                        float(fields.get("Z", "0.0")),
                    ),
                    source=fields.get("SOURCE", "UNKNOWN"),
                )
            )
        except ValueError as exc:
            raise FragmentContractError(f"invalid interaction center row: {row}") from exc
    center_ids = {center.identifier for center in centers}
    interactions: list[AtomCenterInteractionRecord] = []
    for row in _subsection(section, "INTERACTIONS"):
        if row.strip().upper() == "NONE":
            continue
        parts = row.split()
        fields = _key_values(parts[1:])
        center_id = fields.get("CENTER", "")
        if center_id not in center_ids:
            raise FragmentContractError(f"interaction references unknown center: {row}")
        try:
            interactions.append(
                AtomCenterInteractionRecord(
                    identifier=parts[0],
                    kind=fields.get("KIND", "ATOM_CENTER"),
                    atom=int(fields["ATOM"]),
                    center_id=center_id,
                    score=float(fields.get("SCORE", "1.0")),
                    source=fields.get("SOURCE", "UNKNOWN"),
                )
            )
        except (KeyError, ValueError) as exc:
            raise FragmentContractError(f"invalid interaction row: {row}") from exc
    return InteractionCenterDefinition(
        strategy=_section_value(section, "STRATEGY") or "UNKNOWN",
        centers=tuple(centers),
        interactions=tuple(interactions),
    )


def _require_schema(lines: list[str], section_name: str, schema: str) -> None:
    content = section_content(lines, section_name)
    if not content:
        raise FragmentContractError(f"missing #{section_name} section")
    expected = f"SCHEMA {schema}"
    if content[0].strip() != expected:
        raise FragmentContractError(
            f"#{section_name} must start with {expected!r}; found {content[0]!r}"
        )


def _topology_bonds(lines: list[str], *, natoms: int) -> tuple[tuple[int, int], ...]:
    topology = section_content(lines, "TOPOLOGY")
    bond_lines = _subsection(topology, "BONDS")
    bonds: list[tuple[int, int]] = []
    for line in bond_lines:
        if line.strip().upper() == "NONE":
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            i, j = int(parts[0]), int(parts[1])
        except ValueError as exc:
            raise FragmentContractError(f"invalid #TOPOLOGY bond line: {line}") from exc
        if i == j or i < 1 or j < 1 or i > natoms or j > natoms:
            raise FragmentContractError(f"invalid #TOPOLOGY bond indexes: {line}")
        bonds.append(tuple(sorted((i, j))))
    return tuple(sorted(set(bonds)))


def _topology_rings(lines: list[str], *, natoms: int) -> tuple[tuple[int, tuple[int, ...]], ...]:
    topology = section_content(lines, "TOPOLOGY")
    ring_lines = _subsection(topology, "RINGS")
    rings: list[tuple[int, tuple[int, ...]]] = []
    for line in ring_lines:
        if line.strip().upper() == "NONE":
            continue
        parts = line.replace(",", " ").replace("[", " ").replace("]", " ").split()
        if not parts:
            continue
        try:
            ring_index = int(parts[0])
        except ValueError:
            continue
        atoms: list[int] = []
        reading_atoms = False
        for part in parts[1:]:
            token = part.strip()
            if token.upper().startswith("ATOMS="):
                reading_atoms = True
                token = token.split("=", 1)[1]
            elif "=" in token and reading_atoms:
                break
            if not reading_atoms or not token:
                continue
            try:
                atoms.append(int(token))
            except ValueError as exc:
                raise FragmentContractError(f"invalid #TOPOLOGY ring line: {line}") from exc
        if len(atoms) < 3:
            continue
        if any(atom < 1 or atom > natoms for atom in atoms):
            raise FragmentContractError(f"invalid #TOPOLOGY ring atom indexes: {line}")
        rings.append((ring_index, tuple(dict.fromkeys(atoms))))
    return tuple(rings)


def _atom_center_interactions(
    centers: tuple[InteractionCenterRecord, ...],
    *,
    bonds: tuple[tuple[int, int], ...],
    coords: np.ndarray,
    natoms: int,
) -> tuple[AtomCenterInteractionRecord, ...]:
    bonded = {tuple(sorted(bond)) for bond in bonds}
    interactions: list[AtomCenterInteractionRecord] = []
    for center in centers:
        center_atoms = set(center.atoms)
        for atom in range(1, natoms + 1):
            if atom in center_atoms:
                continue
            if any(tuple(sorted((atom, member))) in bonded for member in center_atoms):
                continue
            score = _atom_center_score(atom, center, coords)
            if score <= 0.0:
                continue
            interactions.append(
                AtomCenterInteractionRecord(
                    identifier=f"I{len(interactions) + 1:03d}",
                    kind=f"ATOM_{center.kind}",
                    atom=atom,
                    center_id=center.identifier,
                    score=score,
                    source="AUTO_EQUIDISTANT_GEOMETRY",
                )
            )
    return tuple(interactions)


def _atom_center_score(atom: int, center: InteractionCenterRecord, coords: np.ndarray) -> float:
    atom_coord = coords[atom - 1]
    member_coords = coords[[member - 1 for member in center.atoms]]
    distances = np.linalg.norm(member_coords - atom_coord, axis=1)
    mean_distance = float(np.mean(distances))
    if mean_distance <= RANK_TOLERANCE:
        return 0.0
    spread = float((np.max(distances) - np.min(distances)) / mean_distance)
    center_distance = float(np.linalg.norm(atom_coord - np.asarray(center.center)))
    if center.kind == "RING_CENTER":
        if spread > 0.12 or center_distance > 4.0:
            return 0.0
        return max(0.0, 1.0 - spread / 0.12) * max(0.0, 1.0 - center_distance / 4.0)
    if center.kind == "BOND_CENTER":
        if spread > 0.08 or center_distance > 3.5:
            return 0.0
        return max(0.0, 1.0 - spread / 0.08) * max(0.0, 1.0 - center_distance / 3.5)
    return 0.0


def _connected_components(
    bonds: tuple[tuple[int, int], ...],
    *,
    natoms: int,
) -> tuple[tuple[int, ...], ...]:
    adjacency = {idx: set() for idx in range(1, natoms + 1)}
    for i, j in bonds:
        adjacency[i].add(j)
        adjacency[j].add(i)
    seen: set[int] = set()
    components: list[tuple[int, ...]] = []
    for atom in range(1, natoms + 1):
        if atom in seen:
            continue
        stack = [atom]
        seen.add(atom)
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in sorted(adjacency[current]):
                if neighbor in seen:
                    continue
                seen.add(neighbor)
                stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return tuple(sorted(components, key=lambda item: item[0]))


def _center(coords: np.ndarray, atoms: tuple[int, ...]) -> tuple[float, float, float]:
    center = np.mean(coords[[atom - 1 for atom in atoms]], axis=0)
    return tuple(float(value) for value in center)


def _frame(
    coords: np.ndarray,
    atoms: tuple[int, ...],
) -> tuple[tuple[float, float, float], ...]:
    if _frame_rank(coords, atoms) < 2:
        return _identity_frame()
    p_atom, q_atom = _frame_anchor_atoms(coords, atoms)
    center = np.asarray(_center(coords, atoms), dtype=float)
    p_axis = _unit(coords[p_atom - 1] - center)
    q_axis = _unit(np.cross(p_axis, coords[q_atom - 1] - center))
    s_axis = _unit(np.cross(p_axis, q_axis))
    frame = np.column_stack([p_axis, q_axis, s_axis])
    return tuple(tuple(float(value) for value in frame[:, axis]) for axis in range(3))


def _frame_rank(coords: np.ndarray, atoms: tuple[int, ...]) -> int:
    if len(atoms) < 2:
        return 0
    centered = coords[[atom - 1 for atom in atoms]] - np.asarray(_center(coords, atoms))
    singular_values = np.linalg.svd(centered, compute_uv=False)
    return int(np.sum(singular_values > RANK_TOLERANCE))


def _frame_anchor_atoms(
    coords: np.ndarray,
    atoms: tuple[int, ...],
) -> tuple[int, int]:
    center = np.asarray(_center(coords, atoms), dtype=float)
    ranked = sorted(
        atoms,
        key=lambda atom: (-float(np.linalg.norm(coords[atom - 1] - center)), atom),
    )
    p_atom = ranked[0]
    p_axis = _unit(coords[p_atom - 1] - center)
    q_candidates = []
    for atom in atoms:
        if atom == p_atom:
            continue
        vector = coords[atom - 1] - center
        norm = float(np.linalg.norm(vector))
        if norm <= RANK_TOLERANCE:
            continue
        dot = abs(float(np.dot(p_axis, vector / norm)))
        q_candidates.append((dot, -norm, atom))
    if not q_candidates:
        return ranked[0], ranked[1]
    _dot, _norm, q_atom = min(q_candidates)
    return p_atom, q_atom


def _unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm <= RANK_TOLERANCE:
        raise FragmentContractError("cannot normalize zero-length fragment frame vector")
    return vector / norm


def _identity_frame() -> tuple[tuple[float, float, float], ...]:
    return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


def _subsection(section_lines: list[str], name: str) -> list[str]:
    header = f"[{name.upper()}]"
    start = None
    for idx, line in enumerate(section_lines):
        if line.strip().upper() == header:
            start = idx + 1
            break
    if start is None:
        return []
    end = len(section_lines)
    for idx in range(start, len(section_lines)):
        text = section_lines[idx].strip()
        if text.startswith("[") and text.endswith("]"):
            end = idx
            break
    return list(section_lines[start:end])


def _section_value(section_lines: list[str], key: str) -> str | None:
    key_upper = key.upper()
    for line in section_lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0].upper() == key_upper:
            return parts[1].upper()
    return None


def _key_values(parts: list[str]) -> dict[str, str]:
    fields: dict[str, str] = {}
    for part in parts:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        fields[key.upper()] = value
    return fields


def _parse_int_list(text: str) -> tuple[int, ...]:
    if not text:
        return ()
    try:
        return tuple(int(item) for item in text.replace(";", ",").split(",") if item)
    except ValueError as exc:
        raise FragmentContractError(f"invalid integer list: {text}") from exc


def _center_rows(rows: list[str]) -> dict[str, tuple[float, float, float]]:
    centers: dict[str, tuple[float, float, float]] = {}
    for row in rows:
        parts = row.split()
        if not parts:
            continue
        fields = _key_values(parts[1:])
        try:
            centers[parts[0]] = (
                float(fields.get("X", "0.0")),
                float(fields.get("Y", "0.0")),
                float(fields.get("Z", "0.0")),
            )
        except ValueError as exc:
            raise FragmentContractError(f"invalid center row: {row}") from exc
    return centers


def _frame_rows(rows: list[str]) -> dict[str, tuple[tuple[float, float, float], ...]]:
    frames: dict[str, tuple[tuple[float, float, float], ...]] = {}
    for row in rows:
        parts = row.split()
        if not parts:
            continue
        fields = _key_values(parts[1:])
        axes = []
        for label in ("X", "Y", "Z"):
            text = fields.get(label, "")
            values = text.split(",")
            if len(values) != 3:
                axes = []
                break
            axes.append(tuple(float(value) for value in values))
        if axes:
            frames[parts[0]] = tuple(axes)  # type: ignore[assignment]
    return frames
