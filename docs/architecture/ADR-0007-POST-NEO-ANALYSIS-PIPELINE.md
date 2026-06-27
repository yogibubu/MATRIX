# ADR-0007: Post-NEO Analysis Pipeline

Date: 2026-06-26

## Status

Accepted

## Context

After validation, GICForge constructs the frozen coordinate model for the
molecule. Several workflows can then start:

- GF/PED and vibrational analysis from a Cartesian Hessian;
- rotational analysis and optional SEFit/MORPHEUS;
- thermochemistry;
- anharmonic workflows such as QFF, VPT2 and VCI.

Many of these workflows need data from quantum-mechanical program outputs.
Those outputs must be read through ORACLE adapters, not separately by each
tool.

## Decision

GICForge is the first post-validation producer. It writes frozen `#GIC` and
optional `#SYCART` sections. Downstream tools consume those sections and the
shared QM adapters.

`matrix-gaussian` owns Gaussian file-format I/O:

- Gaussian input writing from enriched XYZ plus `#GIC`;
- Gaussian log/FCHK/QFF/Hessian adapters as they are migrated;
- normalized data models for downstream tools.

The user-facing GICForge workflow may request Gaussian input, but the actual
file serialization is delegated to `matrix-gaussian` to avoid duplicate
Gaussian writers.

Initial command:

```bash
python -m matrix gicforge gaussian-input molecule.xyzin job.gjf \
  --route "#p b3lyp/def2svp opt freq"
```

Downstream ownership:

- `matrix-gf`: consumes frozen `#GIC` plus Cartesian Hessian data.
- `matrix-morpheus`: consumes rotational data, corrections and frozen
  coordinate models for SEFit.
- `matrix-thermo`: consumes normalized frequencies, rotational constants and
  electronic energies.
- `matrix-vpt2-vci`: consumes normalized QFF/anharmonic force-field data.
- `matrix-qm`: owns `#CARTESIAN_HESSIAN`, `#NORMAL_MODES` and `#QFF` as the
  shared handoff sections between QM adapters and scientific tools.

## Consequences

- No downstream package reads Gaussian text directly when an ORACLE adapter
  exists.
- GICForge can produce Gaussian input as a workflow operation without owning a
  private Gaussian writer.
- Cartesian Hessian, normal-mode, QFF, rotational, thermochemical and
  anharmonic data should be promoted to schema-controlled ORACLE sections
  before use by scientific tools.
