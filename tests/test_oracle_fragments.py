from __future__ import annotations

import pytest

from oracle_core import section_content
from oracle_fragments import FragmentContractError, write_fragment_plan_section


def test_fragment_plan_requires_topology_and_synthons(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(FragmentContractError, match="missing #SYNTHONS"):
        write_fragment_plan_section(path)


def test_fragment_plan_section_preserves_existing_sections(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "2",
                "h2",
                "H 0 0 0",
                "H 0 0 1",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
                "INDEXING ATOMS=ONE_BASED",
                "[BONDS]",
                "1 2",
                "",
                "#SYNTHONS",
                "SCHEMA oracle.xyz.synthons.v1",
                "INDEXING ATOMS=ONE_BASED",
                "",
                "#GIC",
                "SCHEMA oracle.xyz.gic.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    write_fragment_plan_section(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    assert section_content(lines, "GIC")[0] == "SCHEMA oracle.xyz.gic.v1"
    fragments = section_content(lines, "FRAGMENTS")
    assert fragments[0] == "SCHEMA oracle.xyz.fragments.v1"
    assert (
        "DEPENDENCIES TOPOLOGY=oracle.xyz.topology.v1 SYNTHONS=oracle.xyz.synthons.v1"
        in fragments
    )
    assert "PENDING ROBUST_TOPOLOGY_CONTRACT" in fragments
