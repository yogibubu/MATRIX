from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class WorkflowStatus(str, Enum):
    MISSING = "missing"
    READY = "ready"
    WARNING = "warning"
    COMPLETE = "complete"


@dataclass(frozen=True)
class WorkflowActionSpec:
    key: str
    label: str
    command: str
    required_sections: tuple[str, ...] = ()
    produced_sections: tuple[str, ...] = ()


@dataclass(frozen=True)
class WindowSpec:
    key: str
    title: str
    description: str
    required_sections: tuple[str, ...] = ()
    produced_sections: tuple[str, ...] = ()
    actions: tuple[WorkflowActionSpec, ...] = ()


ORACLE_GUI_WINDOWS: tuple[WindowSpec, ...] = (
    WindowSpec(
        key="dashboard",
        title="ORACLE Project Dashboard",
        description="Project state, molecule preview, xyzin sections and workflow status.",
        actions=(
            WorkflowActionSpec(
                key="validate",
                label="Validate molecule",
                command="validate",
            ),
        ),
    ),
    WindowSpec(
        key="babel",
        title="ORACLE-Babel / Preprocessing",
        description="Import external geometries and materialize the shared enriched XYZ state.",
        produced_sections=("SOURCE", "BASIC", "SYMMETRY", "TOPOLOGY", "SYNTHONS"),
        actions=(
            WorkflowActionSpec(
                key="preprocess",
                label="Import geometry",
                command="babel preprocess",
            ),
        ),
    ),
    WindowSpec(
        key="avogadro",
        title="Molecule Editor / Avogadro Bridge",
        description="Open the XYZ block in Avogadro and import edited coordinates.",
        required_sections=("SOURCE",),
        actions=(
            WorkflowActionSpec(
                key="replace_xyz",
                label="Import edited XYZ block",
                command="replace xyz block",
            ),
        ),
    ),
    WindowSpec(
        key="topology",
        title="Topology / Synthons / Fragments",
        description="Inspect bonds, rings, charges, synthons and reusable fragments.",
        required_sections=("TOPOLOGY", "SYNTHONS"),
        produced_sections=("FRAGMENTS",),
        actions=(
            WorkflowActionSpec(
                key="fragment_plan",
                label="Plan fragments",
                command="fragments plan",
                required_sections=("TOPOLOGY", "SYNTHONS"),
                produced_sections=("FRAGMENTS",),
            ),
            WorkflowActionSpec(
                key="fragment_build",
                label="Build fragments",
                command="fragments build",
                required_sections=("TOPOLOGY", "SYNTHONS"),
                produced_sections=("FRAGMENTS",),
            ),
        ),
    ),
    WindowSpec(
        key="gicforge",
        title="GICForge",
        description="Build, symmetrize and diagnose frozen GICs and B matrices.",
        required_sections=("SYMMETRY", "TOPOLOGY", "SYNTHONS"),
        produced_sections=("GIC", "SYCART"),
        actions=(
            WorkflowActionSpec(
                key="gic_build",
                label="Build GICs",
                command="gicforge build",
                required_sections=("SYMMETRY", "TOPOLOGY", "SYNTHONS"),
                produced_sections=("GIC",),
            ),
            WorkflowActionSpec(
                key="gic_bmatrix",
                label="Evaluate B matrix",
                command="gicforge bmatrix",
                required_sections=("GIC",),
            ),
            WorkflowActionSpec(
                key="gic_gaussian_input",
                label="Write Gaussian input",
                command="gicforge gaussian-input",
                required_sections=("GIC",),
            ),
        ),
    ),
    WindowSpec(
        key="gf",
        title="GF / PED",
        description="Run harmonic GF/PED from a Cartesian Hessian and frozen GICs.",
        required_sections=("GIC", "CARTESIAN_HESSIAN"),
        produced_sections=("GF_PED",),
        actions=(
            WorkflowActionSpec(
                key="gf_run",
                label="Run GF/PED",
                command="gf",
                required_sections=("GIC", "CARTESIAN_HESSIAN"),
                produced_sections=("GF_PED",),
            ),
        ),
    ),
    WindowSpec(
        key="sefit",
        title="SEFit / MORPHEUS",
        description="Fit semiexperimental structures and multi-molecule MORPHEUS models.",
        required_sections=("ISOTOPOLOGUES",),
        produced_sections=("MORPHEUS",),
        actions=(
            WorkflowActionSpec(
                key="semiexp_fit",
                label="Run SEFit",
                command="semiexp",
                required_sections=("ISOTOPOLOGUES",),
                produced_sections=("MORPHEUS",),
            ),
            WorkflowActionSpec(
                key="semiexp_benchmark",
                label="Run paper benchmark",
                command="semiexp-benchmark",
            ),
        ),
    ),
    WindowSpec(
        key="rovib_thermo",
        title="Thermo / Rovib",
        description="Inspect rotational and vibrational state and run thermochemistry.",
        required_sections=("BASIC", "ROTATIONAL"),
        produced_sections=("THERMO",),
        actions=(
            WorkflowActionSpec(
                key="thermo_run",
                label="Run Thermo",
                command="thermo",
                required_sections=("BASIC", "ROTATIONAL"),
                produced_sections=("THERMO",),
            ),
        ),
    ),
    WindowSpec(
        key="anharmonic",
        title="Anharmonic: VPT2 / VCI / DVR",
        description="Prepare, run and collect VPT2/VCI and DVR workflow state.",
        produced_sections=("VPT2_VCI", "DVR"),
        actions=(
            WorkflowActionSpec(
                key="vpt2_vci_run",
                label="Run VPT2/VCI",
                command="vpt2-vci",
                required_sections=("QFF",),
                produced_sections=("VPT2_VCI",),
            ),
            WorkflowActionSpec(
                key="dvr_run",
                label="Run DVR",
                command="dvr run",
                required_sections=("DVR",),
                produced_sections=("DVR",),
            ),
        ),
    ),
    WindowSpec(
        key="qm_jobs",
        title="QM Jobs",
        description="Generate inputs, monitor external jobs and normalize QM outputs.",
        produced_sections=("CARTESIAN_HESSIAN", "NORMAL_MODES", "QFF"),
        actions=(
            WorkflowActionSpec(
                key="gaussian_from_gic",
                label="Gaussian from GICs",
                command="gicforge gaussian-input",
                required_sections=("GIC",),
            ),
        ),
    ),
    WindowSpec(
        key="diagnostics",
        title="Diagnostics / Regression",
        description="Run numerical audits against corpus, Fortran77 and benchmark fixtures.",
        actions=(
            WorkflowActionSpec(
                key="gic_fortran_audit",
                label="GICForge Python/Fortran audit",
                command="gicforge fortran-audit",
            ),
            WorkflowActionSpec(
                key="gic_corpus_audit",
                label="GIC regression corpus audit",
                command="gicforge corpus-audit",
            ),
        ),
    ),
)


WINDOWS_BY_KEY: dict[str, WindowSpec] = {spec.key: spec for spec in ORACLE_GUI_WINDOWS}


def window_spec(key: str) -> WindowSpec:
    try:
        return WINDOWS_BY_KEY[key]
    except KeyError as exc:
        raise KeyError(f"unknown ORACLE GUI window: {key}") from exc


def all_known_sections() -> tuple[str, ...]:
    sections: set[str] = set()
    for spec in ORACLE_GUI_WINDOWS:
        sections.update(spec.required_sections)
        sections.update(spec.produced_sections)
        for action in spec.actions:
            sections.update(action.required_sections)
            sections.update(action.produced_sections)
    return tuple(sorted(sections))
