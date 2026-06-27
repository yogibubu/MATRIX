from __future__ import annotations

from pathlib import Path
import re

from matrix_chem import read_enriched_xyz
from matrix_core import read_sectioned_lines, section_content


ORACLE_GAUSSIAN_GIC_INPUT_SCHEMA = "oracle.gaussian.gic_input.v1"
REQUIRED_GIC_SCHEMA = "oracle.xyz.gic.v1"
DEFAULT_GIC_ROUTE = "#p hf/sto-3g opt=readallgic"
_GEOM_PAREN_RE = re.compile(r"\bgeom\s*=\s*\((?P<body>[^)]*)\)", flags=re.IGNORECASE)
_GEOM_VALUE_RE = re.compile(
    r"\bgeom\s*=\s*(?P<value>[A-Za-z][A-Za-z0-9_-]*)",
    flags=re.IGNORECASE,
)
_OPT_PAREN_RE = re.compile(r"\bopt\s*=\s*\((?P<body>[^)]*)\)", flags=re.IGNORECASE)
_OPT_VALUE_RE = re.compile(r"\bopt\s*=\s*(?P<value>[A-Za-z][A-Za-z0-9_-]*)", flags=re.IGNORECASE)
_OPT_BARE_RE = re.compile(r"\bopt\b(?!\s*=)", flags=re.IGNORECASE)
_READALLGIC_OPTION_RE = re.compile(r"\breadallgic\b", flags=re.IGNORECASE)
_GICALLSYM_RE = re.compile(r"\bgic(?:all)?symm?\b", flags=re.IGNORECASE)


class GaussianWriteError(ValueError):
    """Raised when ORACLE state cannot be exported to Gaussian input."""


def write_gicforge_gaussian_input(
    enriched_xyz: Path,
    output: Path,
    *,
    route: str = DEFAULT_GIC_ROUTE,
    title: str | None = None,
    charge: int | None = None,
    multiplicity: int | None = None,
    link0: tuple[str, ...] = (),
) -> Path:
    """Write a Gaussian input file from an enriched XYZ carrying #GIC."""
    source = Path(enriched_xyz)
    _require_gic_section(source)
    geometry = read_enriched_xyz(source)
    job_charge = (
        charge if charge is not None else geometry.charge if geometry.charge is not None else 0
    )
    job_multiplicity = (
        multiplicity
        if multiplicity is not None
        else geometry.multiplicity
        if geometry.multiplicity is not None
        else 1
    )
    lines = [
        *[item.strip() for item in link0 if item.strip()],
        _normalize_route(route),
        "",
        title or geometry.comment or source.stem,
        "",
        f"{job_charge} {job_multiplicity}",
    ]
    for atom, (x, y, z) in zip(geometry.atoms, geometry.coordinates_angstrom):
        lines.append(f"{atom:2s} {x:15.8f} {y:15.8f} {z:15.8f}")
    lines.append("")
    gic_lines = _gaussian_gic_lines(source)
    if gic_lines:
        lines.extend(gic_lines)
        lines.append("")
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines), encoding="utf-8")
    return target


def _require_gic_section(path: Path) -> None:
    lines = read_sectioned_lines(Path(path))
    gic = section_content(lines, "GIC")
    if not gic:
        raise GaussianWriteError("missing #GIC section")
    expected = f"SCHEMA {REQUIRED_GIC_SCHEMA}"
    if gic[0].strip() != expected:
        raise GaussianWriteError(f"#GIC must start with {expected!r}; found {gic[0]!r}")


def _gaussian_gic_lines(path: Path) -> list[str]:
    try:
        from matrix_neo import gaussian_gic_lines_from_xyzin
    except ImportError:
        return []
    return gaussian_gic_lines_from_xyzin(Path(path))


def _normalize_route(route: str) -> str:
    text = route.strip()
    if not text:
        raise GaussianWriteError("Gaussian route cannot be empty")
    if not text.startswith("#"):
        text = f"# {text}"
    if _GICALLSYM_RE.search(text):
        raise GaussianWriteError(
            "MATRIX exports already symmetrized GICs; do not request Gaussian GICAllSym"
        )
    if _route_has_opt(text):
        text = _strip_readallgic_geom(text)
        return _collapse_route(_ensure_readallgic_opt(text))
    return _collapse_route(_ensure_readallgic_geom(text))


def _route_has_opt(route: str) -> bool:
    return bool(
        _OPT_PAREN_RE.search(route) or _OPT_VALUE_RE.search(route) or _OPT_BARE_RE.search(route)
    )


def _collapse_route(route: str) -> str:
    return " ".join(route.split())


def _strip_readallgic_geom(route: str) -> str:
    text = _GEOM_PAREN_RE.sub(_remove_readallgic_geom_parenthesized, route)
    return _GEOM_VALUE_RE.sub(_remove_readallgic_geom_value, text)


def _remove_readallgic_geom_parenthesized(match: re.Match[str]) -> str:
    options = [item.strip() for item in match.group("body").split(",") if item.strip()]
    kept = [option for option in options if not _READALLGIC_OPTION_RE.fullmatch(option)]
    if len(kept) == len(options):
        return match.group(0)
    if not kept:
        return ""
    return f"geom=({','.join(kept)})"


def _remove_readallgic_geom_value(match: re.Match[str]) -> str:
    value = match.group("value")
    if value.lower() in {"gic", "readgic", "readallgic"}:
        return ""
    return match.group(0)


def _ensure_readallgic_opt(route: str) -> str:
    if _OPT_PAREN_RE.search(route):
        return _OPT_PAREN_RE.sub(_readallgic_opt_parenthesized, route, count=1)
    if match := _OPT_VALUE_RE.search(route):
        value = match.group("value")
        replacement = (
            "opt=readallgic"
            if value.lower() in {"gic", "readgic", "readallgic"}
            else f"opt=(readallgic,{value})"
        )
        return _OPT_VALUE_RE.sub(replacement, route, count=1)
    if _OPT_BARE_RE.search(route):
        return _OPT_BARE_RE.sub("opt=readallgic", route, count=1)
    return f"{route} opt=readallgic"


def _readallgic_opt_parenthesized(match: re.Match[str]) -> str:
    options = [item.strip() for item in match.group("body").split(",") if item.strip()]
    if not any(_READALLGIC_OPTION_RE.fullmatch(option) for option in options):
        options.insert(0, "readallgic")
    return f"opt=({','.join(options)})"


def _ensure_readallgic_geom(route: str) -> str:
    if _GEOM_PAREN_RE.search(route):
        return _GEOM_PAREN_RE.sub(_readallgic_geom_parenthesized, route, count=1)
    if match := _GEOM_VALUE_RE.search(route):
        value = match.group("value")
        replacement = (
            "geom=readallgic"
            if value.lower() in {"gic", "readgic", "readallgic"}
            else f"geom=(readallgic,{value})"
        )
        return _GEOM_VALUE_RE.sub(replacement, route, count=1)
    return f"{route} geom=readallgic"


def _readallgic_geom_parenthesized(match: re.Match[str]) -> str:
    options = [item.strip() for item in match.group("body").split(",") if item.strip()]
    if not any(_READALLGIC_OPTION_RE.fullmatch(option) for option in options):
        options.insert(0, "readallgic")
    return f"geom=({','.join(options)})"
