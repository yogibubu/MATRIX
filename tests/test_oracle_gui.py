from __future__ import annotations

from oracle_dvr import DVRRequest, dvr_section_from_request, write_dvr_section
from oracle_gui import (
    dvr_gui_state_lines,
    load_dvr_gui_state,
    load_vpt2_vci_gui_state,
    vpt2_vci_gui_state_lines,
)
from oracle_vpt2_vci import vpt2_vci_section_from_run, write_vpt2_vci_section


def test_gui_dvr_state_reads_refreshed_dvr_section(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    outdir = tmp_path / "dvr"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=outdir,
        figdir=tmp_path / "fig",
        prefix="demo",
    )
    write_dvr_section(xyzin, dvr_section_from_request(request))
    outdir.mkdir()
    (outdir / "demo_summary.txt").write_text("Mass-weighted path Hamiltonian\n", encoding="utf-8")
    (outdir / "demo_levels.csv").write_text(
        "state,energy_cm-1,energy_above_ground_cm-1\n0,25.0,0.0\n",
        encoding="utf-8",
    )
    (outdir / "demo_grid.csv").write_text(
        "grid,s_au,s_sqrtamu_angstrom,V_cm-1\n0,0.0,0.0,0.0\n",
        encoding="utf-8",
    )

    state = load_dvr_gui_state(xyzin, refresh=True)
    lines = dvr_gui_state_lines(state)

    assert state.ready
    assert state.level_count == 1
    assert "status: complete" in lines
    assert "ready: 1" in lines


def test_gui_vpt2_vci_state_reads_refreshed_section(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    run_dir = tmp_path / "vpt2_vci"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")
    run_dir.mkdir()
    (run_dir / "vpt2_vci.report").write_text("report\n", encoding="utf-8")
    (run_dir / "vpt2_vci_frequencies.csv").write_text(
        "mode,harmonic_frequency_cm-1\n1,100.0\n",
        encoding="utf-8",
    )
    (run_dir / "vpt2_vci_comparison.csv").write_text(
        "root,vpt2_abs_cm-1,vci_abs_cm-1,delta_abs_cm-1,"
        "vpt2_exc_cm-1,vci_exc_cm-1,delta_exc_cm-1\n"
        "1,50.0,50.0,0.0,0.0,0.0,0.0\n",
        encoding="utf-8",
    )
    (run_dir / "vpt2_vci_mode_contributions.csv").write_text(
        "root,mode,expected_quanta\n1,1,0.5\n",
        encoding="utf-8",
    )
    write_vpt2_vci_section(
        xyzin,
        vpt2_vci_section_from_run(
            source_kind="xyzin",
            source_path=xyzin,
            run_dir=run_dir,
            report_path=run_dir / "vpt2_vci.report",
            csv_dir=run_dir,
            manifest_path=None,
            max_quanta=2,
            roots=1,
            vci_method="dense",
            status="prepared",
        ),
    )

    state = load_vpt2_vci_gui_state(xyzin, refresh=True)
    lines = vpt2_vci_gui_state_lines(state)

    assert state.ready
    assert state.root_count == 1
    assert state.mode_count == 1
    assert "status: complete" in lines
    assert "ready: 1" in lines
