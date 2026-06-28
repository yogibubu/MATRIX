from __future__ import annotations

from pathlib import Path

from matrix_core.cli import main
from matrix_molpro import molpro_job_status, run_molpro_job
from matrix_oracle.qm_jobs import OracleQMJobsController
from matrix_orca import orca_job_status, run_orca_job


def test_molpro_launcher_uses_shared_external_qm_runner(tmp_path):
    (tmp_path / "molpro.com").write_text("***,water\ngeometry={H}\n", encoding="utf-8")
    executable = _fake_executable(
        tmp_path / "fake_molpro.sh",
        "printf 'Molpro calculation terminated\\nVariable memory released\\n'\n",
    )

    result = run_molpro_job(tmp_path, executable=str(executable))
    status = molpro_job_status(tmp_path)

    assert result.success is True
    assert result.exit_code == 0
    assert result.output_path == tmp_path / "molpro.out"
    assert status.status == "completed"
    assert status.normal_termination is True
    assert "Variable memory released" in result.output_path.read_text(encoding="utf-8")


def test_orca_launcher_uses_shared_external_qm_runner(tmp_path):
    (tmp_path / "orca.inp").write_text("! HF STO-3G\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
    executable = _fake_executable(
        tmp_path / "fake_orca.sh",
        "printf '****ORCA TERMINATED NORMALLY****\\n'\n",
    )

    result = run_orca_job(tmp_path, executable=str(executable))
    status = orca_job_status(tmp_path)

    assert result.success is True
    assert result.exit_code == 0
    assert result.output_path == tmp_path / "orca.out"
    assert status.status == "completed"
    assert status.normal_termination is True


def test_molpro_and_orca_cli_run_commands(tmp_path, capsys):
    molpro_dir = tmp_path / "molpro"
    molpro_dir.mkdir()
    (molpro_dir / "molpro.com").write_text("***,h\ngeometry={H}\n", encoding="utf-8")
    molpro_exe = _fake_executable(
        tmp_path / "fake_molpro.sh",
        "printf 'Molpro calculation terminated\\n'\n",
    )

    assert main(["molpro", "run", str(molpro_dir), "--executable", str(molpro_exe)]) == 0
    molpro_output = capsys.readouterr().out
    assert "program: Molpro" in molpro_output
    assert "success: 1" in molpro_output

    orca_dir = tmp_path / "orca"
    orca_dir.mkdir()
    (orca_dir / "orca.inp").write_text("! SP\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
    orca_exe = _fake_executable(
        tmp_path / "fake_orca.sh",
        "printf 'ORCA TERMINATED NORMALLY\\n'\n",
    )

    assert main(["orca", "run", str(orca_dir), "--executable", str(orca_exe)]) == 0
    orca_output = capsys.readouterr().out
    assert "program: ORCA" in orca_output
    assert "success: 1" in orca_output


def test_qm_jobs_controller_exposes_molpro_and_orca_launchers(tmp_path):
    controller = OracleQMJobsController(tmp_path / "molecule.xyzin")

    molpro_run = controller.molpro_run_command(
        tmp_path / "molpro",
        executable="molpro",
        input_path="input.com",
        background=True,
    )
    orca_status = controller.orca_status_command(tmp_path / "orca", output_path="orca.out")

    assert molpro_run.argv[1:4] == ("-m", "matrix", "molpro")
    assert "run" in molpro_run.argv
    assert "--background" in molpro_run.argv
    assert orca_status.argv[1:4] == ("-m", "matrix", "orca")
    assert "status" in orca_status.argv
    assert "--output" in orca_status.argv


def _fake_executable(path: Path, body: str) -> Path:
    path.write_text(
        "#!/bin/sh\n"
        'test -f "$1" || exit 9\n'
        f"{body}"
        "exit 0\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path
