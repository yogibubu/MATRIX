from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SectionCompletionHint:
    section: str
    window: str
    action: str
    command: str
    note: str = ""

    def line(self) -> str:
        text = f"#{self.section}: use {self.window} -> {self.action} ({self.command})"
        if self.note:
            text += f"; {self.note}"
        return text


SECTION_COMPLETION_HINTS: dict[str, SectionCompletionHint] = {
    "SOURCE": SectionCompletionHint(
        "SOURCE",
        "Structure / ORACLE-Babel",
        "Preprocess",
        "oracle babel preprocess",
    ),
    "BASIC": SectionCompletionHint(
        "BASIC",
        "Structure / ORACLE-Babel",
        "Preprocess",
        "oracle babel preprocess",
    ),
    "SYMMETRY": SectionCompletionHint(
        "SYMMETRY",
        "Structure / ORACLE-Babel",
        "Preprocess",
        "oracle babel preprocess",
    ),
    "TOPOLOGY": SectionCompletionHint(
        "TOPOLOGY",
        "Structure / ORACLE-Babel",
        "Preprocess",
        "oracle babel preprocess",
        "topology is computed once from the normalized geometry",
    ),
    "SYNTHONS": SectionCompletionHint(
        "SYNTHONS",
        "Structure / ORACLE-Babel",
        "Preprocess",
        "oracle babel preprocess",
        "synthons are derived from the shared topology/charge model",
    ),
    "FRAGMENTS": SectionCompletionHint(
        "FRAGMENTS",
        "Structure",
        "Build Fragments",
        "oracle fragments build",
    ),
    "GIC": SectionCompletionHint(
        "GIC",
        "GICForge",
        "Build GICs",
        "oracle gicforge build",
    ),
    "SYCART": SectionCompletionHint(
        "SYCART",
        "GICForge",
        "Build GICs with SYCART",
        "oracle gicforge build --sycart",
    ),
    "CARTESIAN_HESSIAN": SectionCompletionHint(
        "CARTESIAN_HESSIAN",
        "QM Jobs",
        "Promote Gaussian FCHK",
        "oracle gaussian promote-fchk",
        "the GF/PED tab can also accept an FCHK file directly",
    ),
    "NORMAL_MODES": SectionCompletionHint(
        "NORMAL_MODES",
        "QM Jobs",
        "Promote Gaussian FCHK",
        "oracle gaussian promote-fchk",
    ),
    "QFF": SectionCompletionHint(
        "QFF",
        "QM Jobs",
        "Promote Gaussian FCHK with QFF data",
        "oracle gaussian promote-fchk",
    ),
    "ROTATIONAL": SectionCompletionHint(
        "ROTATIONAL",
        "QM Jobs",
        "Promote rovibrational Gaussian log",
        "oracle gaussian promote-rovib",
    ),
    "VIBRATIONAL": SectionCompletionHint(
        "VIBRATIONAL",
        "QM Jobs",
        "Promote rovibrational Gaussian log",
        "oracle gaussian promote-rovib",
    ),
    "DELTABVIB": SectionCompletionHint(
        "DELTABVIB",
        "QM Jobs",
        "Promote rovibrational Gaussian log",
        "oracle gaussian promote-rovib",
    ),
    "GF_PED": SectionCompletionHint(
        "GF_PED",
        "GF/PED",
        "Run GF/PED",
        "oracle gf",
    ),
    "ISOTOPOLOGUES": SectionCompletionHint(
        "ISOTOPOLOGUES",
        "SEFit / Rotational Spectroscopy",
        "Import or define isotopologues",
        "oracle semiexp --job ... --xyzin ...",
        "the #ISOTOPOLOGUES section is the standalone SEFit input contract",
    ),
    "MORPHEUS": SectionCompletionHint(
        "MORPHEUS",
        "SEFit",
        "Run SEFit",
        "oracle semiexp",
    ),
    "VPT2_VCI": SectionCompletionHint(
        "VPT2_VCI",
        "Anharmonic",
        "Run or collect VPT2/VCI",
        "oracle vpt2-vci",
    ),
    "DVR": SectionCompletionHint(
        "DVR",
        "Anharmonic",
        "Prepare, run or collect DVR",
        "oracle dvr prepare/run/collect",
    ),
    "THERMO": SectionCompletionHint(
        "THERMO",
        "Thermo/Kinetics",
        "Run Thermo",
        "oracle thermo",
    ),
    "KINETICS": SectionCompletionHint(
        "KINETICS",
        "Thermo/Kinetics",
        "Prepare kinetics workflow",
        "oracle kinetics",
        "planned section; use normalized thermo/rovib outputs as inputs",
    ),
    "TRINITY": SectionCompletionHint(
        "TRINITY",
        "TRINITY",
        "Prepare TRINITY",
        "oracle trinity prepare",
    ),
}


def section_completion_hint(section: str) -> SectionCompletionHint | None:
    return SECTION_COMPLETION_HINTS.get(_normalize_section(section))


def missing_sections_guidance(sections: tuple[str, ...] | list[str]) -> tuple[str, ...]:
    lines: list[str] = []
    for section in sections:
        hint = section_completion_hint(section)
        if hint is None:
            lines.append(f"#{_normalize_section(section)}: add this section to xyzin with the owning adapter/tool")
        else:
            lines.append(hint.line())
    return tuple(lines)


def missing_sections_message(sections: tuple[str, ...] | list[str]) -> str:
    normalized = tuple(_normalize_section(section) for section in sections)
    if not normalized:
        return "project is not available"
    lines = ["missing " + ", ".join(normalized), "Suggested completion:"]
    lines.extend(f"- {line}" for line in missing_sections_guidance(normalized))
    return "\n".join(lines)


def _normalize_section(section: str) -> str:
    return section.strip().upper().lstrip("#")
