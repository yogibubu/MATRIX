from __future__ import annotations

from io import StringIO

from matrix_chem import AtomicSynthons, build_topology_objects
from matrix_chem.topology.topology_writer import write_topology_section
from matrix_neo.definition import _topology_bonds, _topology_rings, topology_bond_orders_from_lines


def test_topology_builds_graph_and_atomic_synthons_for_water():
    coords = [
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 1.0, 0.0),
    ]
    atomic_numbers = [8, 1, 1]

    continuous, discrete, ringset, synthons, aromaticity = build_topology_objects(
        coords,
        atomic_numbers,
    )

    assert continuous.natoms == 3
    assert sorted(discrete.bonds) == [(0, 1), (0, 2)]
    assert ringset.rings == []
    assert isinstance(synthons, AtomicSynthons)
    assert synthons.Zeff(0) > 0.0
    assert synthons.canonical_signature(0)[0] == 8
    assert aromaticity.aromatic_atoms == set()


def test_topology_writer_uses_modern_xyzin_contract_for_neo_consumers():
    coords = [
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 1.0, 0.0),
    ]
    atomic_numbers = [8, 1, 1]
    continuous, discrete, ringset, synthons, aromaticity = build_topology_objects(
        coords,
        atomic_numbers,
    )

    stream = StringIO()
    write_topology_section(
        stream,
        cg=continuous,
        dg=discrete,
        ringset=ringset,
        synthons=synthons,
        aromaticity=aromaticity,
    )

    lines = [line.rstrip("\n") for line in stream.getvalue().splitlines()]
    topology_start = lines.index("#TOPOLOGY") + 1
    topology = lines[topology_start:]
    assert topology[0] == "SCHEMA oracle.xyz.topology.v1"
    assert "[BONDS]" in topology
    assert "[BOND_ORDERS]" in topology
    assert "[RINGS]" in topology
    assert _topology_bonds(lines, natoms=3) == ((1, 2), (1, 3))
    assert topology_bond_orders_from_lines(lines, natoms=3).keys() == {(1, 2), (1, 3)}
    assert _topology_rings(lines, natoms=3) == ()
