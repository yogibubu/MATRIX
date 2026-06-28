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


MOLPRO_EXECUTABLE = "molpro"
MOLPRO_SPEC = ExternalQMProgramSpec(
    name="Molpro",
    default_executable=MOLPRO_EXECUTABLE,
    env_var="MATRIX_MOLPRO_EXE",
    input_candidates=("molpro.com", "molpro.inp", "input.com", "input.inp"),
    output_suffix=".out",
    normal_markers=("Molpro calculation terminated", "Variable memory released"),
    error_markers=("ERROR EXIT", "Molpro calculation terminated with errors"),
    pid_file="molpro.pid",
)


def molpro_job_status(
    workdir: Path | str,
    *,
    input_path: Path | str | None = None,
    output_path: Path | str | None = None,
) -> ExternalQMJobStatus:
    return external_qm_job_status(
        workdir,
        MOLPRO_SPEC,
        input_path=input_path,
        output_path=output_path,
    )


def run_molpro_job(
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
        MOLPRO_SPEC,
        executable=executable,
        input_path=input_path,
        output_path=output_path,
        background=background,
        timeout=timeout,
        env=env,
        extra_args=extra_args,
    )
