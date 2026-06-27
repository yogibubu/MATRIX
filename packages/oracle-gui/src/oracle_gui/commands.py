from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shlex
import subprocess
import sys
from collections.abc import Sequence


@dataclass(frozen=True)
class OracleGuiCommand:
    label: str
    argv: tuple[str, ...]
    required_sections: tuple[str, ...] = ()
    produced_sections: tuple[str, ...] = ()
    cwd: Path | None = None

    def shell_line(self) -> str:
        return " ".join(shlex.quote(item) for item in self.argv)

    def run(self, *, timeout: float | None = None) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            self.argv,
            cwd=None if self.cwd is None else str(self.cwd),
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )


def validate_command(xyzin: Path | str, *, require_fragments: bool = False) -> OracleGuiCommand:
    argv = [*_oracle(), "validate", str(Path(xyzin))]
    if require_fragments:
        argv.append("--require-fragments")
    return OracleGuiCommand("Validate molecule", tuple(argv))


def avogadro_command(
    xyzin: Path | str,
    *,
    executable: str = "avogadro2",
) -> OracleGuiCommand:
    return OracleGuiCommand("Open in Avogadro", (executable, str(Path(xyzin))))


def preprocess_command(
    source: Path | str,
    output: Path | str,
    *,
    source_kind: str = "auto",
    symmetry_distance: float | None = None,
    symmetry_inertia: float | None = None,
    max_rotation_order: int | None = None,
) -> OracleGuiCommand:
    argv = [*_oracle(), "babel", "preprocess", str(Path(source)), str(Path(output))]
    if source_kind != "auto":
        argv.extend(["--source-kind", source_kind])
    if symmetry_distance is not None:
        argv.extend(["--symmetry-distance", str(symmetry_distance)])
    if symmetry_inertia is not None:
        argv.extend(["--symmetry-inertia", str(symmetry_inertia)])
    if max_rotation_order is not None:
        argv.extend(["--max-rotation-order", str(max_rotation_order)])
    return OracleGuiCommand(
        "Preprocess molecule",
        tuple(argv),
        produced_sections=("SOURCE", "BASIC", "SYMMETRY", "TOPOLOGY", "SYNTHONS"),
    )


def fragments_command(xyzin: Path | str, action: str = "build") -> OracleGuiCommand:
    if action not in {"plan", "build", "centers"}:
        raise ValueError(f"unsupported fragments action: {action}")
    return OracleGuiCommand(
        f"Fragments {action}",
        (*_oracle(), "fragments", action, str(Path(xyzin))),
        required_sections=("TOPOLOGY", "SYNTHONS"),
        produced_sections=("FRAGMENTS",),
    )


def gicforge_build_command(
    xyzin: Path | str,
    *,
    symmetrize: bool = True,
    sycart: bool = True,
    improper_dihedrals: bool = True,
) -> OracleGuiCommand:
    argv = [*_oracle(), "gicforge", "build", str(Path(xyzin))]
    _append_flag(argv, "--symmetrize", symmetrize)
    _append_flag(argv, "--sycart", sycart)
    _append_flag(argv, "--improper-dihedrals", improper_dihedrals)
    return OracleGuiCommand(
        "Build GICForge coordinates",
        tuple(argv),
        required_sections=("SYMMETRY", "TOPOLOGY", "SYNTHONS"),
        produced_sections=("GIC", "SYCART"),
    )


def gicforge_report_command(xyzin: Path | str, output: Path | str | None = None) -> OracleGuiCommand:
    argv = [*_oracle(), "gicforge", "report", str(Path(xyzin))]
    if output is not None:
        argv.append(str(Path(output)))
    return OracleGuiCommand("Write GICForge report", tuple(argv), required_sections=("GIC",))


def gicforge_bmatrix_command(xyzin: Path | str, output: Path | str | None = None) -> OracleGuiCommand:
    argv = [*_oracle(), "gicforge", "bmatrix", str(Path(xyzin))]
    if output is not None:
        argv.append(str(Path(output)))
    return OracleGuiCommand("Evaluate GIC B matrix", tuple(argv), required_sections=("GIC",))


def gicforge_gaussian_input_command(
    xyzin: Path | str,
    output: Path | str,
    *,
    route: str = "#p hf/sto-3g",
    title: str | None = None,
    charge: int | None = None,
    multiplicity: int | None = None,
) -> OracleGuiCommand:
    argv = [*_oracle(), "gicforge", "gaussian-input", str(Path(xyzin)), str(Path(output))]
    argv.extend(["--route", route])
    if title is not None:
        argv.extend(["--title", title])
    if charge is not None:
        argv.extend(["--charge", str(charge)])
    if multiplicity is not None:
        argv.extend(["--multiplicity", str(multiplicity)])
    return OracleGuiCommand("Write Gaussian GIC input", tuple(argv), required_sections=("GIC",))


def gf_command(
    xyzin: Path | str,
    *,
    fchk: Path | str | None = None,
    out: Path | str | None = None,
    csv_dir: Path | str | None = None,
    local: bool = False,
    symmetry_blocks: bool = True,
    force_threshold: float | None = None,
) -> OracleGuiCommand:
    argv = [*_oracle(), "gf", "--xyzin", str(Path(xyzin))]
    if fchk is not None:
        argv.extend(["--fchk", str(Path(fchk))])
    if out is not None:
        argv.extend(["--out", str(Path(out))])
    if csv_dir is not None:
        argv.extend(["--csv-dir", str(Path(csv_dir))])
    _append_flag(argv, "--local", local)
    _append_flag(argv, "--symmetry-blocks", symmetry_blocks)
    if force_threshold is not None:
        argv.extend(["--force-threshold", str(force_threshold)])
    return OracleGuiCommand(
        "Run GF/PED",
        tuple(argv),
        required_sections=("GIC",),
        produced_sections=("GF_PED",),
    )


def thermo_command(xyzin: Path | str, *, out: Path | str | None = None) -> OracleGuiCommand:
    argv = [*_oracle(), "thermo", str(Path(xyzin))]
    if out is not None:
        argv.extend(["--out", str(Path(out))])
    return OracleGuiCommand(
        "Run Thermo",
        tuple(argv),
        required_sections=("BASIC", "ROTATIONAL"),
        produced_sections=("THERMO",),
    )


def vpt2_vci_command(
    xyzin: Path | str,
    *,
    run_dir: Path | str | None = None,
    max_quanta: int = 2,
    roots: int = 10,
) -> OracleGuiCommand:
    argv = [
        *_oracle(),
        "vpt2-vci",
        "--xyzin",
        str(Path(xyzin)),
        "--max-quanta",
        str(max_quanta),
        "--roots",
        str(roots),
    ]
    if run_dir is not None:
        argv.extend(["--run-dir", str(Path(run_dir))])
    return OracleGuiCommand(
        "Run VPT2/VCI",
        tuple(argv),
        required_sections=("QFF",),
        produced_sections=("VPT2_VCI",),
    )


def vpt2_vci_collect_command(xyzin: Path | str, *, no_write: bool = False) -> OracleGuiCommand:
    argv = [*_oracle(), "vpt2-vci", "--collect", str(Path(xyzin))]
    _append_flag(argv, "--no-write", no_write)
    return OracleGuiCommand(
        "Collect VPT2/VCI outputs",
        tuple(argv),
        produced_sections=("VPT2_VCI",),
    )


def dvr_run_command(
    xyzin: Path | str,
    *,
    timeout: float | None = None,
    check_only: bool = False,
) -> OracleGuiCommand:
    argv = [*_oracle(), "dvr", "run", "--xyzin", str(Path(xyzin))]
    _append_flag(argv, "--check-only", check_only)
    if timeout is not None:
        argv.extend(["--timeout", str(timeout)])
    return OracleGuiCommand("Run DVR", tuple(argv), required_sections=("DVR",), produced_sections=("DVR",))


def dvr_collect_command(xyzin: Path | str, *, no_write: bool = False) -> OracleGuiCommand:
    argv = [*_oracle(), "dvr", "collect", str(Path(xyzin))]
    _append_flag(argv, "--no-write", no_write)
    return OracleGuiCommand("Collect DVR outputs", tuple(argv), produced_sections=("DVR",))


def semiexp_command(
    job: Path | str,
    outdir: Path | str,
    *,
    backend: str = "python",
    extra_args: Sequence[str] = (),
) -> OracleGuiCommand:
    argv = [*_oracle(), "semiexp", "--job", str(Path(job)), "--outdir", str(Path(outdir)), "--backend", backend]
    argv.extend(str(item) for item in extra_args)
    return OracleGuiCommand(
        "Run SEFit / MORPHEUS",
        tuple(argv),
        required_sections=("ISOTOPOLOGUES",),
        produced_sections=("MORPHEUS",),
    )


def semiexp_benchmark_command(
    *,
    snapshot: Path | str | None = None,
    outdir: Path | str | None = None,
    refresh: bool = True,
    update_snapshot: bool = False,
) -> OracleGuiCommand:
    argv = [*_oracle(), "semiexp-benchmark"]
    if snapshot is not None:
        argv.extend(["--snapshot", str(Path(snapshot))])
    if outdir is not None:
        argv.extend(["--outdir", str(Path(outdir))])
    _append_flag(argv, "--no-refresh", not refresh)
    _append_flag(argv, "--update-snapshot", update_snapshot)
    return OracleGuiCommand("Run SEFit/MORPHEUS paper benchmark", tuple(argv))


def gicforge_fortran_audit_command(
    *,
    root: Path | str | None = None,
    workdir: Path | str | None = None,
    molecules: Sequence[str] = (),
    limit: int | None = None,
    tolerance: float = 2.0e-8,
    output_format: str = "summary",
) -> OracleGuiCommand:
    argv = [*_oracle(), "gicforge", "fortran-audit"]
    if root is not None:
        argv.extend(["--root", str(Path(root))])
    if workdir is not None:
        argv.extend(["--workdir", str(Path(workdir))])
    for molecule in molecules:
        argv.extend(["--molecule", molecule])
    if limit is not None:
        argv.extend(["--limit", str(limit)])
    argv.extend(["--tolerance", str(tolerance), "--format", output_format])
    return OracleGuiCommand("Run GICForge Python/Fortran audit", tuple(argv))


def _oracle() -> tuple[str, str, str]:
    return (sys.executable, "-m", "oracle")


def _append_flag(argv: list[str], flag: str, enabled: bool) -> None:
    if enabled:
        argv.append(flag)
