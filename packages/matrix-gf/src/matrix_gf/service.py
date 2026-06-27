from __future__ import annotations

from dataclasses import dataclass
import csv
from io import StringIO
from pathlib import Path
import re

import numpy as np

from matrix_neo import read_gic_definition_from_xyzin
from matrix_gaussian import hessian_input_from_gaussian_fchk
from matrix_qm import hessian_input_from_xyzin

from .internal import (
    GFLocalOptions,
    InternalGFResult,
    gf_from_hessian_input_and_xyzin,
    gf_from_hessian_input_with_matrix_gics,
    topology_bonds_from_xyzin,
)
from .nonbonded import nonbonded_cartesian_hessian_correction, synthon_charges_from_xyzin


@dataclass(frozen=True)
class GFReport:
    fchk_path: Path
    result: InternalGFResult
    text: str
    xyzin_path: Path | None = None
    scale_path: Path | None = None
    hessian_source: str = ""


@dataclass(frozen=True)
class GFScalingClass:
    """Named class of GICs receiving the same diagonal Pulay factor."""

    name: str
    factor: float
    patterns: tuple[str, ...]

    def validate(self) -> None:
        if not self.name.strip():
            raise ValueError("Pulay scaling class name cannot be empty")
        if not self.patterns or any(not pattern.strip() for pattern in self.patterns):
            raise ValueError(f"Pulay scaling class {self.name!r} needs at least one pattern")
        if self.factor < 0.0:
            raise ValueError("Pulay scaling factors must be non-negative")


@dataclass(frozen=True)
class GFScalingRulePreview:
    order: int
    kind: str
    name: str
    factor: float
    patterns: tuple[str, ...]
    matches: tuple[int, ...]
    family: str = ""


@dataclass(frozen=True)
class GFScalingAssignment:
    index: int
    identifier: str
    name: str
    family: str
    factor: float
    source: str
    label: str


@dataclass(frozen=True)
class GFScalingPreview:
    xyzin_path: Path | None
    assignments: tuple[GFScalingAssignment, ...]
    rules: tuple[GFScalingRulePreview, ...]

    @property
    def changed_count(self) -> int:
        return sum(1 for item in self.assignments if abs(item.factor - 1.0) > 1.0e-12)


@dataclass(frozen=True)
class _ScalingEvent:
    kind: str
    selector: str = ""
    factor: float = 1.0
    scaling_class: GFScalingClass | None = None


def run_gf_report_from_fchk(fchk_path: Path) -> GFReport:
    """Read an FCHK adapter and return a formatted quick GF/PED report."""
    path = Path(fchk_path)
    hessian_input = hessian_input_from_gaussian_fchk(path)
    result = gf_from_hessian_input_with_matrix_gics(hessian_input)
    return GFReport(path, result, format_gf_report(path, result))


def run_xyzin_gf_report_from_fchk(
    fchk_path: Path,
    xyzin_path: Path,
    *,
    scale_path: Path | None = None,
    scale_records: tuple[str, ...] = (),
    scale_class_records: tuple[str, ...] = (),
    local: bool = False,
    force_threshold: float | None = None,
    block_by_irrep: bool = False,
    subtract_electrostatic: bool = False,
    subtract_uff_vdw: bool = False,
    nonbonded_14_scale: float = 0.5,
) -> GFReport:
    """Run the frozen-xyzin GF branch from a Cartesian Hessian FCHK adapter."""
    path = Path(fchk_path)
    hessian_input = hessian_input_from_gaussian_fchk(path)
    return _run_xyzin_gf_report_from_hessian_input(
        hessian_input,
        path,
        Path(xyzin_path),
        hessian_source=f"FCHK {path}",
        scale_path=scale_path,
        scale_records=scale_records,
        scale_class_records=scale_class_records,
        local=local,
        force_threshold=force_threshold,
        block_by_irrep=block_by_irrep,
        subtract_electrostatic=subtract_electrostatic,
        subtract_uff_vdw=subtract_uff_vdw,
        nonbonded_14_scale=nonbonded_14_scale,
    )


def run_xyzin_gf_report_from_xyzin(
    xyzin_path: Path,
    *,
    scale_path: Path | None = None,
    scale_records: tuple[str, ...] = (),
    scale_class_records: tuple[str, ...] = (),
    local: bool = False,
    force_threshold: float | None = None,
    block_by_irrep: bool = False,
    subtract_electrostatic: bool = False,
    subtract_uff_vdw: bool = False,
    nonbonded_14_scale: float = 0.5,
) -> GFReport:
    """Run GF from frozen #GIC and #CARTESIAN_HESSIAN sections in one xyzin."""
    xyzin = Path(xyzin_path)
    hessian_input = hessian_input_from_xyzin(xyzin)
    return _run_xyzin_gf_report_from_hessian_input(
        hessian_input,
        xyzin,
        xyzin,
        hessian_source=f"#CARTESIAN_HESSIAN in {xyzin}",
        scale_path=scale_path,
        scale_records=scale_records,
        scale_class_records=scale_class_records,
        local=local,
        force_threshold=force_threshold,
        block_by_irrep=block_by_irrep,
        subtract_electrostatic=subtract_electrostatic,
        subtract_uff_vdw=subtract_uff_vdw,
        nonbonded_14_scale=nonbonded_14_scale,
    )


def _run_xyzin_gf_report_from_hessian_input(
    hessian_input,
    source_path: Path,
    xyzin: Path,
    *,
    hessian_source: str,
    scale_path: Path | None = None,
    scale_records: tuple[str, ...] = (),
    scale_class_records: tuple[str, ...] = (),
    local: bool = False,
    force_threshold: float | None = None,
    block_by_irrep: bool = False,
    subtract_electrostatic: bool = False,
    subtract_uff_vdw: bool = False,
    nonbonded_14_scale: float = 0.5,
) -> GFReport:
    xyzin = Path(xyzin)
    definition = read_gic_definition_from_xyzin(xyzin)
    correction = None
    correction_label = "NONE"
    if subtract_electrostatic or subtract_uff_vdw:
        topology_bonds = topology_bonds_from_xyzin(xyzin)
        charges = None
        labels: list[str] = []
        if subtract_electrostatic:
            charges, charge_source = synthon_charges_from_xyzin(xyzin, len(hessian_input.atomic_numbers))
            labels.append(f"ELECTROSTATIC({charge_source})")
        if subtract_uff_vdw:
            labels.append("UFF_VDW")
        correction = nonbonded_cartesian_hessian_correction(
            hessian_input.cartesian_coordinates_bohr,
            hessian_input.atomic_numbers,
            topology_bonds,
            charges=charges,
            electrostatic=subtract_electrostatic,
            uff_vdw=subtract_uff_vdw,
            one_four_scale=nonbonded_14_scale,
        )
        correction_label = "+".join(labels) + f"; 1-4 scale={float(nonbonded_14_scale):g}"
    names = tuple(gic.name for gic in definition.gics)
    labels = tuple(_gic_display_label(gic.identifier, gic.gaussian_expression) for gic in definition.gics)
    scaling = pulay_scaling_factors(
        len(definition.gics),
        labels=labels,
        names=names,
        scale_path=scale_path,
        scale_records=scale_records,
        scale_class_records=scale_class_records,
    )
    result = gf_from_hessian_input_and_xyzin(
        hessian_input,
        xyzin,
        scaling_factors=scaling,
        local_options=GFLocalOptions(enabled=local),
        force_threshold=force_threshold,
        block_by_irrep=block_by_irrep,
        cartesian_hessian_correction=correction,
        cartesian_hessian_correction_label=correction_label,
    )
    return GFReport(
        Path(source_path),
        result,
        format_gf_report(
            Path(source_path),
            result,
            xyzin_path=xyzin,
            scale_path=scale_path,
            hessian_source=hessian_source,
        ),
        xyzin_path=xyzin,
        scale_path=scale_path,
        hessian_source=hessian_source,
    )


def format_gf_report(
    fchk_path: Path,
    result: InternalGFResult,
    *,
    xyzin_path: Path | None = None,
    scale_path: Path | None = None,
    hessian_source: str | None = None,
) -> str:
    lines = [
        "GF/PED from ORACLE non-redundant GICs",
        f"Hessian source: {hessian_source or f'FCHK {Path(fchk_path)}'}",
        f"Coordinate source: {result.coordinate_source}",
        f"Point group: {result.point_group}",
        f"Symmetrized GICs: {result.symmetrized_gics}",
        f"Matrix model: {result.matrix_model}",
        f"Cartesian Hessian correction: {result.hessian_correction}",
        f"GIC count: {len(result.gic_labels)}",
        "",
        "Frequencies (cm-1):",
    ]
    if xyzin_path is not None:
        lines.insert(2, f"Frozen xyzin: {Path(xyzin_path)}")
    if result.scaling_factors is not None:
        changed = int(np.sum(np.abs(result.scaling_factors - 1.0) > 1.0e-12))
        lines.insert(
            3 if xyzin_path is not None else 2,
            "Pulay Hessian scaling: "
            f"applied ({changed} factors != 1; file={Path(scale_path) if scale_path else 'inline/default'})",
        )
    if result.force_threshold is not None:
        lines.insert(
            3 if xyzin_path is not None else 2,
            f"Force-constant threshold: {result.force_threshold:g}",
        )
    if result.block_labels:
        counts = ", ".join(
            f"{label}:{result.block_labels.count(label)}"
            for label in dict.fromkeys(result.block_labels)
        )
        lines.insert(3 if xyzin_path is not None else 2, f"GF blocks: {counts}")
    for idx, freq in enumerate(result.frequencies_cm, start=1):
        lines.append(f"  mode {idx:3d}: {freq:12.3f}")

    lines.extend(["", "GIC labels:"])
    for idx, label in enumerate(result.gic_labels, start=1):
        name = result.gic_names[idx - 1] if idx <= len(result.gic_names) else f"GIC{idx:03d}"
        irrep = result.gic_irreps[idx - 1] if idx <= len(result.gic_irreps) else "UNK"
        lines.append(f"  GIC{idx:03d}: {name:12s} irrep={irrep:6s} {label}")

    lines.extend(["", "PED (%) rows=GIC cols=modes:"])
    header = "          " + " ".join(f"M{idx:02d}" for idx in range(1, len(result.frequencies_cm) + 1))
    lines.append(header)
    for idx, row in enumerate(result.ped.values, start=1):
        values = " ".join(f"{value:7.2f}" for value in row)
        lines.append(f"  GIC{idx:03d} {values}")
    return "\n".join(lines)


def gf_csv_tables(report: GFReport) -> dict[str, str]:
    """Return CSV tables for GF frequencies, GIC labels, matrices, normal modes and PED."""
    freq_rows = [["mode", "frequency_cm-1"]]
    freq_rows.extend(
        [[idx, f"{freq:.10g}"] for idx, freq in enumerate(report.result.frequencies_cm, start=1)]
    )

    label_rows = [["gic", "name", "irrep", "label"]]
    for idx, label in enumerate(report.result.gic_labels, start=1):
        name = report.result.gic_names[idx - 1] if idx <= len(report.result.gic_names) else f"GIC{idx:03d}"
        irrep = report.result.gic_irreps[idx - 1] if idx <= len(report.result.gic_irreps) else "UNK"
        label_rows.append([f"GIC{idx:03d}", name, irrep, label])

    ped_rows = [
        ["gic", "name", "irrep", *[f"mode_{idx}" for idx in range(1, len(report.result.frequencies_cm) + 1)]]
    ]
    for idx, row in enumerate(report.result.ped.values, start=1):
        name = report.result.gic_names[idx - 1] if idx <= len(report.result.gic_names) else f"GIC{idx:03d}"
        irrep = report.result.gic_irreps[idx - 1] if idx <= len(report.result.gic_irreps) else "UNK"
        ped_rows.append([f"GIC{idx:03d}", name, irrep, *[f"{value:.10g}" for value in row]])

    return {
        "frequencies.csv": _csv_text(freq_rows),
        "gic_labels.csv": _csv_text(label_rows),
        "ped.csv": _csv_text(ped_rows),
        "normal_modes.csv": _csv_text(_gic_mode_table(report.result.modes_internal, "mode", report.result)),
        "force_constants.csv": _csv_text(_square_gic_table(report.result.force_constants, report.result)),
        "g_matrix.csv": _csv_text(_square_gic_table(report.result.g_matrix, report.result)),
    }


def write_csv_tables(report: GFReport, outdir: Path, *, prefix: str = "gf") -> dict[str, Path]:
    target_dir = Path(outdir)
    target_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for name, text in gf_csv_tables(report).items():
        path = target_dir / f"{prefix}_{name}"
        path.write_text(text, encoding="utf-8")
        written[name] = path
    return written


def pulay_scaling_factors(
    n_gics: int,
    *,
    labels: tuple[str, ...],
    names: tuple[str, ...] = (),
    scale_path: Path | None = None,
    scale_records: tuple[str, ...] = (),
    scale_class_records: tuple[str, ...] = (),
) -> np.ndarray | None:
    """Build diagonal Pulay scaling factors from selector and class records."""
    events = _collect_scaling_events(
        scale_path=scale_path,
        scale_records=scale_records,
        scale_class_records=scale_class_records,
    )
    if not events:
        return None
    factors = np.ones(n_gics, dtype=float)
    for event in events:
        if event.kind == "selector":
            if event.factor < 0.0:
                raise ValueError("Pulay scaling factors must be non-negative")
            for index in _resolve_scaling_selector(event.selector, labels=labels, names=names):
                factors[index] = event.factor
            continue
        if event.kind == "class" and event.scaling_class is not None:
            event.scaling_class.validate()
            for index in _resolve_scaling_class(event.scaling_class, labels=labels, names=names):
                factors[index] = event.scaling_class.factor
            continue
        raise ValueError(f"Unsupported Pulay scaling event: {event!r}")
    return factors


def gf_scaling_preview_from_xyzin(
    xyzin_path: Path,
    *,
    scale_path: Path | None = None,
    scale_records: tuple[str, ...] = (),
    scale_class_records: tuple[str, ...] = (),
) -> GFScalingPreview:
    """Return a dry-run Pulay scaling assignment preview for frozen xyzin GICs."""
    xyzin = Path(xyzin_path)
    definition = read_gic_definition_from_xyzin(xyzin)
    names = tuple(gic.name for gic in definition.gics)
    labels = tuple(_gic_display_label(gic.identifier, gic.gaussian_expression) for gic in definition.gics)
    return pulay_scaling_preview(
        len(definition.gics),
        labels=labels,
        names=names,
        xyzin_path=xyzin,
        scale_path=scale_path,
        scale_records=scale_records,
        scale_class_records=scale_class_records,
    )


def pulay_scaling_preview(
    n_gics: int,
    *,
    labels: tuple[str, ...],
    names: tuple[str, ...] = (),
    xyzin_path: Path | None = None,
    scale_path: Path | None = None,
    scale_records: tuple[str, ...] = (),
    scale_class_records: tuple[str, ...] = (),
) -> GFScalingPreview:
    """Resolve Pulay scaling records without running GF."""
    events = _collect_scaling_events(
        scale_path=scale_path,
        scale_records=scale_records,
        scale_class_records=scale_class_records,
    )
    factors = np.ones(n_gics, dtype=float)
    sources = ["default" for _ in range(n_gics)]
    padded_names = _padded_names(names, n_gics)
    rule_previews: list[GFScalingRulePreview] = []
    for order, event in enumerate(events, start=1):
        if event.kind == "selector":
            if event.factor < 0.0:
                raise ValueError("Pulay scaling factors must be non-negative")
            matches = _resolve_scaling_selector(event.selector, labels=labels, names=names)
            source = f"selector {event.selector}"
            for index in matches:
                factors[index] = event.factor
                sources[index] = source
            rule_previews.append(
                GFScalingRulePreview(
                    order=order,
                    kind="selector",
                    name=event.selector,
                    factor=event.factor,
                    patterns=(event.selector,),
                    matches=tuple(index + 1 for index in matches),
                    family=_common_family(matches, labels=labels, names=padded_names),
                )
            )
            continue
        if event.kind == "class" and event.scaling_class is not None:
            scaling_class = event.scaling_class
            scaling_class.validate()
            matches = _resolve_scaling_class(scaling_class, labels=labels, names=names)
            source = f"class {scaling_class.name}"
            for index in matches:
                factors[index] = scaling_class.factor
                sources[index] = source
            rule_previews.append(
                GFScalingRulePreview(
                    order=order,
                    kind="class",
                    name=scaling_class.name,
                    factor=scaling_class.factor,
                    patterns=scaling_class.patterns,
                    matches=tuple(index + 1 for index in matches),
                    family=_common_family(matches, labels=labels, names=padded_names),
                )
            )
            continue
        raise ValueError(f"Unsupported Pulay scaling event: {event!r}")
    assignments = tuple(
        GFScalingAssignment(
            index=index + 1,
            identifier=f"GIC{index + 1:03d}",
            name=padded_names[index] or f"GIC{index + 1:03d}",
            family=_gic_coordinate_family(padded_names[index], labels[index]) or "unknown",
            factor=float(factors[index]),
            source=sources[index],
            label=labels[index],
        )
        for index in range(n_gics)
    )
    return GFScalingPreview(None if xyzin_path is None else Path(xyzin_path), assignments, tuple(rule_previews))


def format_gf_scaling_preview(preview: GFScalingPreview) -> str:
    """Format a Pulay scaling dry-run preview for CLI/GUI logs."""
    lines = ["GF/PED Pulay scaling preview"]
    if preview.xyzin_path is not None:
        lines.append(f"Frozen xyzin: {preview.xyzin_path}")
    lines.append(f"GIC count: {len(preview.assignments)}")
    lines.append(f"Changed factors: {preview.changed_count}")
    lines.append("")
    lines.append("Rules:")
    if not preview.rules:
        lines.append("  none; all factors remain 1.0")
    for rule in preview.rules:
        pattern_text = "|".join(rule.patterns)
        family = rule.family or "mixed/unknown"
        matches = ",".join(f"GIC{index:03d}" for index in rule.matches)
        lines.append(
            f"  {rule.order:02d} {rule.kind:8s} {rule.name:20s} "
            f"factor={rule.factor:.8g} family={family:10s} "
            f"matches={len(rule.matches):3d} [{matches}] patterns={pattern_text}"
        )
    lines.append("")
    lines.append("Assignments:")
    lines.append("  GIC     name             family     factor       source               label")
    for item in preview.assignments:
        lines.append(
            f"  {item.identifier:7s} {item.name[:14]:14s} {item.family[:10]:10s} "
            f"{item.factor:11.8g} {item.source[:20]:20s} {item.label}"
        )
    return "\n".join(lines)


def _collect_scaling_events(
    *,
    scale_path: Path | None = None,
    scale_records: tuple[str, ...] = (),
    scale_class_records: tuple[str, ...] = (),
) -> list[_ScalingEvent]:
    events: list[_ScalingEvent] = []
    if scale_path is not None:
        events.extend(_read_scaling_events(Path(scale_path)))
    for record in scale_records:
        events.append(_selector_event(*_parse_scaling_record(record)))
    for record in scale_class_records:
        events.append(_class_event(_parse_scaling_class_record(record)))
    return events


def _read_scaling_records(path: Path) -> list[tuple[str, float]]:
    records: list[tuple[str, float]] = []
    for event in _read_scaling_events(path):
        if event.kind == "selector":
            records.append((event.selector, event.factor))
    return records


def _read_scaling_events(path: Path) -> list[_ScalingEvent]:
    events: list[_ScalingEvent] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.split("#", 1)[0].split("!", 1)[0].strip()
        if not line:
            continue
        lower = line.lower().replace(",", " ").split()
        if lower[:2] in (["selector", "factor"], ["gic", "factor"], ["name", "factor"]):
            continue
        events.append(_parse_scaling_event(line))
    return events


def _parse_scaling_event(record: str) -> _ScalingEvent:
    text = record.strip()
    if text.lower().startswith("class "):
        return _class_event(_parse_scaling_class_record(text))
    return _selector_event(*_parse_scaling_record(text))


def _selector_event(selector: str, factor: float) -> _ScalingEvent:
    return _ScalingEvent("selector", selector=selector, factor=factor)


def _class_event(scaling_class: GFScalingClass) -> _ScalingEvent:
    return _ScalingEvent("class", scaling_class=scaling_class)


def _parse_scaling_record(record: str) -> tuple[str, float]:
    text = record.strip()
    if "=" in text:
        selector, value = text.split("=", 1)
    elif "," in text:
        selector, value = text.split(",", 1)
    else:
        parts = text.split()
        if len(parts) != 2:
            raise ValueError(f"Invalid Pulay scaling record: {record!r}")
        selector, value = parts
    return selector.strip(), float(value.strip())


def _parse_scaling_class_record(record: str) -> GFScalingClass:
    text = record.strip()
    if text.lower().startswith("class "):
        text = text.split(None, 1)[1].strip()
    if ":" in text:
        parts = text.split(":", 2)
        if len(parts) != 3:
            raise ValueError(f"Invalid Pulay scaling class record: {record!r}")
        name, factor_text, patterns_text = parts
    else:
        parts = text.split(None, 2)
        if len(parts) != 3:
            raise ValueError(f"Invalid Pulay scaling class record: {record!r}")
        name, factor_text, patterns_text = parts
    factor_text = factor_text.strip()
    if "=" in factor_text:
        _key, factor_text = factor_text.split("=", 1)
    patterns = tuple(pattern.strip() for pattern in patterns_text.split("|") if pattern.strip())
    scaling_class = GFScalingClass(name=name.strip(), factor=float(factor_text.strip()), patterns=patterns)
    scaling_class.validate()
    return scaling_class


def _resolve_scaling_selector(
    selector: str,
    *,
    labels: tuple[str, ...],
    names: tuple[str, ...],
) -> tuple[int, ...]:
    token = selector.strip()
    lower = token.lower()
    n_gics = len(labels)
    if lower in {"*", "all", "default"}:
        return tuple(range(n_gics))
    if lower.startswith("gic") and lower[3:].isdigit():
        index = int(lower[3:]) - 1
        if 0 <= index < n_gics:
            return (index,)
    if token.isdigit():
        index = int(token) - 1
        if 0 <= index < n_gics:
            return (index,)
    exact = [idx for idx, name in enumerate(names) if name == token]
    exact.extend(idx for idx, label in enumerate(labels) if label == token)
    if exact:
        return tuple(dict.fromkeys(exact))
    matches = [
        idx
        for idx, (name, label) in enumerate(zip(_padded_names(names, n_gics), labels))
        if token in name or token in label
    ]
    if len(matches) == 1:
        return (matches[0],)
    if matches:
        raise ValueError(f"Pulay scaling selector {token!r} is ambiguous")
    raise ValueError(f"Pulay scaling selector {token!r} did not match any GIC")


def _resolve_scaling_class(
    scaling_class: GFScalingClass,
    *,
    labels: tuple[str, ...],
    names: tuple[str, ...],
) -> tuple[int, ...]:
    matches: list[int] = []
    for pattern in scaling_class.patterns:
        matches.extend(_resolve_scaling_class_pattern(pattern, labels=labels, names=names))
    unique = tuple(dict.fromkeys(matches))
    if not unique:
        raise ValueError(f"Pulay scaling class {scaling_class.name!r} did not match any GIC")
    padded_names = _padded_names(names, len(labels))
    families = {
        _gic_coordinate_family(padded_names[index], labels[index])
        for index in unique
    }
    known_families = {family for family in families if family}
    if len(known_families) > 1:
        raise ValueError(
            f"Pulay scaling class {scaling_class.name!r} mixes coordinate types: "
            f"{', '.join(sorted(known_families))}"
        )
    return unique


def _resolve_scaling_class_pattern(
    pattern: str,
    *,
    labels: tuple[str, ...],
    names: tuple[str, ...],
) -> tuple[int, ...]:
    token = pattern.strip()
    lower = token.lower()
    n_gics = len(labels)
    if lower in {"*", "all", "default"}:
        return tuple(range(n_gics))
    if lower.startswith("gic") and lower[3:].isdigit():
        index = int(lower[3:]) - 1
        if 0 <= index < n_gics:
            return (index,)
    if token.isdigit():
        index = int(token) - 1
        if 0 <= index < n_gics:
            return (index,)
    padded_names = _padded_names(names, n_gics)
    exact = [idx for idx, name in enumerate(padded_names) if name.lower() == lower]
    exact.extend(idx for idx, label in enumerate(labels) if label.lower() == lower)
    if exact:
        return tuple(dict.fromkeys(exact))
    matches = [
        idx
        for idx, (name, label) in enumerate(zip(padded_names, labels))
        if lower in name.lower() or lower in label.lower()
    ]
    if matches:
        return tuple(dict.fromkeys(matches))
    raise ValueError(f"Pulay scaling class pattern {token!r} did not match any GIC")


def _gic_coordinate_family(name: str, label: str) -> str:
    text = f"{name} {label}".lower()
    if any(token in text for token in ("rpck", "qpck", "phip", "pck", "butterfly")):
        return "ring"
    if any(token in text for token in ("fragment", "frag", "centroid", "center", "centre")):
        return "special"
    if re.search(r"\br\s*\(", text) or "str" in text or "stretch" in text or "bond(" in text:
        return "stretch"
    if re.search(r"\ba\s*\(", text) or "bend" in text or "angle(" in text:
        return "bend"
    if re.search(r"\bd\s*\(", text) or "tors" in text or "dih" in text or "dihedral(" in text:
        return "torsion"
    if re.search(r"\bu\s*\(", text) or "oop" in text or "improper" in text or "out_of_plane(" in text:
        return "oop"
    if re.search(r"\bl\s*\(", text) or "linear_bend" in text or "lin" in text:
        return "linear"
    return ""


def _common_family(
    indices: tuple[int, ...],
    *,
    labels: tuple[str, ...],
    names: tuple[str, ...],
) -> str:
    families = {
        _gic_coordinate_family(names[index], labels[index])
        for index in indices
    }
    known = {family for family in families if family}
    if len(known) == 1:
        return next(iter(known))
    if len(known) > 1:
        return "mixed"
    return ""


def _padded_names(names: tuple[str, ...], n_gics: int) -> tuple[str, ...]:
    if len(names) >= n_gics:
        return names[:n_gics]
    return names + tuple("" for _ in range(n_gics - len(names)))


def _gic_display_label(identifier: str, expression: str) -> str:
    return f"{identifier} {expression}" if expression and expression != "NONE" else identifier


def _csv_text(rows: list[list[object]]) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerows(rows)
    return stream.getvalue()


def _gic_mode_table(matrix: np.ndarray, label: str, result: InternalGFResult) -> list[list[object]]:
    values = np.asarray(matrix, dtype=float)
    rows: list[list[object]] = [
        ["gic", "name", "irrep", *[f"{label}_{idx}" for idx in range(1, values.shape[1] + 1)]]
    ]
    for idx, row in enumerate(values, start=1):
        name = result.gic_names[idx - 1] if idx <= len(result.gic_names) else f"GIC{idx:03d}"
        irrep = result.gic_irreps[idx - 1] if idx <= len(result.gic_irreps) else "UNK"
        rows.append([f"GIC{idx:03d}", name, irrep, *[f"{value:.10g}" for value in row]])
    return rows


def _square_gic_table(matrix: np.ndarray, result: InternalGFResult) -> list[list[object]]:
    values = np.asarray(matrix, dtype=float)
    rows: list[list[object]] = [
        ["gic", "name", "irrep", *[f"GIC{idx:03d}" for idx in range(1, values.shape[1] + 1)]]
    ]
    for idx, row in enumerate(values, start=1):
        name = result.gic_names[idx - 1] if idx <= len(result.gic_names) else f"GIC{idx:03d}"
        irrep = result.gic_irreps[idx - 1] if idx <= len(result.gic_irreps) else "UNK"
        rows.append([f"GIC{idx:03d}", name, irrep, *[f"{value:.10g}" for value in row]])
    return rows
