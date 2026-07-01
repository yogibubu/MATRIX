# MATRIX Internal Coordinate Framework

This document is the MATRIX adaptation of the internal-coordinate architecture
notes from the `newcoord` repository.  It is a refactoring contract, not a
request to redesign the mathematics.  Merlino/GICForge numerical behaviour is
preserved unless a change is explicitly justified, tested, and documented.

## Guiding Rules

- Topology is built before geometry-dependent analysis.
- Geometry is evaluated before coordinate values and B rows.
- Chemically meaningful coordinate definitions are generated before numerical
  reduction.
- Numerical linear algebra is confined to reduction, symmetrization and GF-like
  stages.
- Every stage appends information to the shared molecular/xyzin model and must
  not silently invalidate earlier stages.
- Deterministic output is mandatory whenever chemically possible.
- New coordinate types are integrated by adding a generator, registering it and
  adding regression tests.

## Pipeline

The canonical coordinate workflow is:

```text
Input
  -> topology perception
  -> topology analysis and validation
  -> geometry evaluation
  -> coordinate generation
  -> coordinate classification
  -> redundancy elimination
  -> symmetrization
  -> validation and serialization
```

In MATRIX this maps to:

| Stage | MATRIX owner | Notes |
| --- | --- | --- |
| Input/import | LINK and QM adapters | Reads XYZ, QM output or SMILES; writes shared xyzin sections. |
| Topology | LINK / matrix-chem | Owns bonds, rings, fragments, synthons and validation. |
| Geometry | LINK / matrix-chem | Owns coordinates, inertia, symmetry thresholds and auxiliary-node positions. |
| Coordinate generation | NEO/GICForge | Builds primitive and chemically adapted coordinates. |
| Classification | NEO policy layer | Assigns coordinate family, reduction class and symmetry block. |
| Redundancy elimination | NEO reducer | Uses analytic B rows and type-local rank logic. |
| Symmetrization | NEO symmetrizer | Projects only after a globally non-redundant basis exists. |
| Downstream use | GF, MORPHEUS, DVR, VPT2/VCI | Consume frozen GIC/SYCART definitions; they do not rebuild coordinates. |

## Data Ownership

Topology contains graph-only data:

- atoms and atom properties;
- bonds and bond orders;
- rings, fused rings, spiro and bridged systems;
- fragments and pseudo-bonds;
- auxiliary-node definitions;
- synthon/equivalence classes.

The topology contract is serialized in `#TOPOLOGY` with schema
`oracle.xyz.topology.v1`.  NEO/GICForge consumes `[BONDS]`, `[BOND_ORDERS]` and
`[RINGS]` from this section; it may report how the data were used, but it must
not rediscover a second graph.  Bond orders are normalized before NEO sees
them: Gaussian Mayer values are preferred when present, otherwise the shared
continuous Pauling model is used.  Ring-puckering flexibility is therefore a
relative intra-ring policy applied to normalized bond orders, not a private
simple/double/aromatic classifier inside the coordinate generator.

Pseudo-bonds are topology policy, not a hidden GIC patch.  In
`SPECIAL_COORDINATES` mode fragments remain separate and NEO adds protected
inter-fragment coordinates.  In `PSEUDO_BONDS` mode selected weak contacts are
inserted into the graph before primitive generation; they can generate ordinary
stretches, bends and torsions, but they must not create artificial ring
coordinates.

Geometry contains coordinate-dependent data:

- Cartesian coordinates;
- centers of mass and geometric centers;
- inertia tensor and principal axes;
- fragment orientations;
- instantiated auxiliary-node positions.

Coordinate definitions contain:

- coordinate type and family;
- participating atoms or auxiliary nodes;
- parent primitives for composite coordinates;
- generating algorithm;
- reduction class;
- symmetry block;
- provenance.

The frozen `#GIC` and `#SYCART` xyzin sections are the contract consumed by
later tools.  GF, MORPHEUS and anharmonic modules must not reparse Gaussian
ReadAllGIC blocks or regenerate topology to reinterpret those coordinates.

## Coordinate Generator Boundary

Each coordinate type has one logical generator.  The current production path is
still the Merlino-compatible NEO/GICForge implementation; the generator registry
is a migration boundary for gradual refactoring.

A generator may:

- read topology and geometry;
- create coordinate definitions;
- attach metadata and provenance.

A generator must not:

- mutate topology;
- perform redundancy elimination;
- perform symmetrization;
- diagonalize matrices or solve eigenproblems;
- manipulate Hessians.

Current extraction status:

| Generator | Status | Behavioural contract |
| --- | --- | --- |
| `StretchGenerator` | Extracted in `matrix_neo.generators.generate_stretch_coordinates` | Must remain identical to the stretch part of `definition._primitive_candidates`. |
| `LocalXHStretchGenerator` | Extracted in `matrix_neo.generators.generate_stretch_coordinates` | Uses the same opt-in X-H policy as NEO; GF treats only these rows as local X-H. |
| `LocalSymmetryAngleSALCGenerator` | Refactor boundary around the Merlino local stretch/angle path | Bond lengths are grouped by endpoint-equivalence and length before local SVD/SALC; valence angles 2-4 follow the Merlino templates guided by local ligand classes; high coordination 5-9 uses recognized polyhedron templates, with ligand-equivalence plus local SVD/SALC as fallback. |
| Other generators | Still embedded in the Merlino-compatible NEO path | Must be extracted one family at a time with golden tests before behaviour changes. |

## GIC Construction Closure Criteria

The GIC construction is considered complete only when these conditions hold for
both Python and Fortran paths where both exist:

- the same normalized topology inputs produce the same primitive families,
  protected coordinate groups, rank and active non-redundant coordinate count;
- ring coordinates are never eliminated by unrelated coordinate families;
- cyclic bends, RPck components, butterfly coordinates, out-of-plane/improper
  coordinates, linear bends, fragment-center coordinates and local X--H
  stretches have explicit reduction classes;
- symmetrization builds SALCs without mixing coordinate families and preserves
  protected special-coordinate classes;
- non-total-symmetric coordinates are frozen in Gaussian optimization inputs,
  while total-symmetric coordinates remain active;
- diagnostics list protected coordinates, eliminated coordinates, rank
  decisions, symmetry blocks and ring bond-order flexibility factors;
- the golden corpus compares MATRIX Python, MATRIX Fortran and original
  Merlino behaviour for labels/families, symmetry blocks and selected B rows.

## Reduction and Symmetry Boundary

Reduction receives generated coordinates and analytic B rows.  It may use
modified Gram-Schmidt, rank-revealing logic, BBt diagnostics or transformed
Hessian information when available.  Tie-breaking must prefer chemically
meaningful coordinates when numerically stable.

Symmetrization receives the already non-redundant coordinate set.  It must
preserve coordinate families and avoid mixing unrelated physical motions.  A
symmetry-adapted coordinate basis that cannot reproduce the expected
vibrational representation is a contract failure, not a downstream GF problem.

## Regression Corpus

Changes to coordinate construction, reduction, symmetry labels, B rows or
downstream GF use require tests on representative systems:

- acyclic molecules;
- isolated rings;
- fused and polycyclic rings;
- spiro and bridged systems;
- linear molecules;
- weak complexes and pseudo-bonds;
- metal/auxiliary-center cases when relevant;
- Python/Fortran comparison where the Fortran path exists.

Golden tests should compare at least coordinate counts, labels/families,
selected primitive/GIC rows, symmetry diagnostics and B-row numerical values
where available.
