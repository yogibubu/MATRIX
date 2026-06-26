from __future__ import annotations

from pathlib import Path

from tools import oracle_run


def test_lcb25_fetch_cli_calls_sync(tmp_path, monkeypatch, capsys):
    calls = {}

    def fake_sync(root, *, datasets=None, force=False):
        calls["root"] = root
        calls["datasets"] = datasets
        calls["force"] = force
        return Path(root) / "manifest.json"

    monkeypatch.setattr("oracle_babel.sync_lcb25_library", fake_sync)

    rc = oracle_run.main(
        ["lcb25", "fetch", "--root", str(tmp_path / "cache"), "--dataset", "se", "--force"]
    )

    assert rc == 0
    assert calls == {"root": tmp_path / "cache", "datasets": ["se"], "force": True}
    assert "Synced LCB25 library" in capsys.readouterr().out


def test_fragments_plan_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyz"

    def fake_write(target):
        calls["target"] = target

    monkeypatch.setattr("oracle_fragments.write_fragment_plan_section", fake_write)

    rc = oracle_run.main(["fragments", "plan", str(path)])

    assert rc == 0
    assert calls == {"target": path}
    assert "Planned ORACLE fragment workflow" in capsys.readouterr().out
