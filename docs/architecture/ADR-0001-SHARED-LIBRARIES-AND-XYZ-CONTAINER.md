# ADR 0001: Shared Libraries And Enriched XYZ Container

Date: 2026-06-26

## Decision

All ORACLE tools must use shared libraries for shared tasks. A module must not
reimplement parsing, topology perception, isotope handling, GIC construction,
Gaussian parsing, manifest writing, section management or backend execution when
an ORACLE library already owns that function.

The main communication object between tools is an enriched XYZ container:

1. The file starts with a valid ordinary XYZ block.
2. Additional information is appended in named uppercase sections.
3. Each tool owns only the section or sections it produces.
4. A tool that updates one section must preserve all unrelated sections exactly.
5. External formats such as Gaussian logs, TOML jobs, CSV tables and MSR inputs
   are adapters or import sources. They must be materialized into the enriched
   XYZ container before downstream ORACLE tools consume them.

This keeps tools independently developable while sharing the same molecule,
metadata and workflow state.

## Rationale

ORACLE is a pipeline, not a collection of disconnected scripts. The same
operations must have one implementation:

- atom and isotope data;
- XYZ and enriched-section parsing;
- molecular graph, topology, ring and symmetry perception;
- GIC and B-matrix construction;
- Gaussian input/output parsing;
- rotational/vibrational correction conventions;
- semiexperimental observation records;
- backend discovery and execution;
- run manifests and checksums.

When every module reads and writes the same enriched XYZ container through the
same library APIs, GICForge, MORPHEUS, GF/PED, DVR, VPT2/VCI and GUI workflows
can be developed independently without duplicating data models.

## Container Shape

Minimum file:

```text
3
water example
O   0.000000   0.000000   0.000000
H   0.000000   0.757000   0.586000
H   0.000000  -0.757000   0.586000
```

Enriched file:

```text
3
water example
O   0.000000   0.000000   0.000000
H   0.000000   0.757000   0.586000
H   0.000000  -0.757000   0.586000

#BASIC
SCHEMA oracle.xyz.basic.v1
CHARGE 0
MULTIPLICITY 1

#TOPOLOGY
SCHEMA oracle.xyz.topology.v1

#GIC
SCHEMA oracle.xyz.gic.v1

#ISOTOPOLOGUES
SCHEMA oracle.xyz.isotopologues.v1
```

## Section Ownership

Initial section ownership:

- `#BASIC`: `oracle-core` and GUI/project importers.
- `#SMILES`: SMILES/import adapters.
- `#GAUSSIAN`: `oracle-gaussian`.
- `#GAUSSIAN_TOPOLOGY`: `oracle-gaussian` as adapter data, consumed by
  `oracle-chem`.
- `#TOPOLOGY`: `oracle-chem`.
- `#SYMMETRY`: `oracle-chem`.
- `#GIC`: `oracle-gicforge`.
- `#SYCART`: `oracle-gicforge`.
- `#ROTATIONAL`: `oracle-chem` for derived constants, `oracle-morpheus` for
  fitted/equilibrium records.
- `#VIBRATIONAL`: `oracle-gaussian`, `oracle-gf` or `oracle-vpt2-vci` depending
  on source, with schema-specific ownership.
- `#ISOTOPOLOGUES`: shared schema owned by `oracle-core`/`oracle-morpheus`.
- `#MORPHEUS`: `oracle-morpheus`.
- `#GF_PED`: `oracle-gf`.
- `#DVR`: `oracle-dvr`.
- `#VPT2_VCI`: `oracle-vpt2-vci`.

When a new section is needed, its schema and owner must be documented before
implementation.

## Library Rule

Each package can expose domain services, but shared primitives live in lower
packages:

```text
oracle-core
  sectioned XYZ file utilities, manifests, workspace, config, errors
oracle-chem
  atom/mass/geometry/topology/ring/symmetry model
oracle-gaussian
  Gaussian adapters only
oracle-gicforge
  GIC/SYCART schemas and B matrices
oracle-morpheus, oracle-gf, oracle-dvr, oracle-vpt2-vci
  consume shared models and append their own sections
oracle-gui
  calls services; does not own scientific logic
```

No package may reach upward for shared logic. For example, `oracle-gicforge`
may use `oracle-chem`, but `oracle-chem` must not import `oracle-gicforge`.

## Consequences

- Every workflow can accept a single enriched XYZ file as canonical input.
- Tools can add value incrementally: import geometry, add topology, add GICs,
  add isotopologues, add MORPHEUS results, add GF/PED tables, and so on.
- Tests must check section preservation whenever a tool updates the container.
- Compatibility readers can accept old `merlino.xyzin.*` schemas, but new
  output should use `oracle.xyz.*` schemas.

