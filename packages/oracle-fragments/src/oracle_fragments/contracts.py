from __future__ import annotations

from pathlib import Path

from oracle_core import read_sectioned_lines, replace_section, section_content


ORACLE_XYZ_FRAGMENTS_SCHEMA = "oracle.xyz.fragments.v1"
ORACLE_XYZ_FRAGMENT_LIBRARY_SCHEMA = "oracle.xyz.fragment_library.v1"
ORACLE_XYZ_ASSEMBLY_SCHEMA = "oracle.xyz.assembly.v1"

REQUIRED_TOPOLOGY_SCHEMA = "oracle.xyz.topology.v1"
REQUIRED_SYNTHONS_SCHEMA = "oracle.xyz.synthons.v1"


class FragmentContractError(ValueError):
    """Raised when a file cannot enter the ORACLE fragment workflow."""


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


def _require_schema(lines: list[str], section_name: str, schema: str) -> None:
    content = section_content(lines, section_name)
    if not content:
        raise FragmentContractError(f"missing #{section_name} section")
    expected = f"SCHEMA {schema}"
    if content[0].strip() != expected:
        raise FragmentContractError(
            f"#{section_name} must start with {expected!r}; found {content[0]!r}"
        )
