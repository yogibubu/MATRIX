from __future__ import annotations

from pathlib import Path
from collections.abc import Mapping, Sequence

from matrix_engines import (
    ExternalQMJobStatus,
    ExternalQMProgramSpec,
    ExternalQMRunResult,
    external_qm_job_status,
    run_external_qm_job,
)


ORCA_EXECUTABLE = "orca"
ORCA_SPEC = ExternalQMProgramSpec(
    name="ORCA",
    default_executable=ORCA_EXECUTABLE,
    env_var="MATRIX_ORCA_EXE",
    input_candidates=("orca.inp", "input.inp", "molecule.inp"),
    output_suffix=".out",
    normal_markers=("ORCA TERMINATED NORMALLY",),
    error_markers=("ORCA finished by error termination", "ORCA TERMINATED ABNORMALLY"),
    pid_file="orca.pid",
)


def orca_job_status(
    workdir: Path | str,
    *,
    input_path: Path | str | None = None,
    output_path: Path | str | None = None,
) -> ExternalQMJobStatus:
    return external_qm_job_status(
        workdir,
        ORCA_SPEC,
        input_path=input_path,
        output_path=output_path,
    )


def run_orca_job(
    workdir: Path | str,
    *,
    executable: str | None = None,
    input_path: Path | str | None = None,
    output_path: Path | str | None = None,
    background: bool = False,
    timeout: float | None = None,
    env: Mapping[str, str] | None = None,
    extra_args: Sequence[str] = (),
) -> ExternalQMRunResult:
    return run_external_qm_job(
        workdir,
        ORCA_SPEC,
        executable=executable,
        input_path=input_path,
        output_path=output_path,
        background=background,
        timeout=timeout,
        env=env,
        extra_args=extra_args,
    )
