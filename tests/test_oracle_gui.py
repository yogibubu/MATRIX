from __future__ import annotations

from oracle_dvr import DVRRequest, dvr_section_from_request, write_dvr_section
from oracle_gui import dvr_gui_state_lines, load_dvr_gui_state


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
