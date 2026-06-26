# ORACLE

**Operational Recognition of Atomistic Connectivity and Local Environments**

ORACLE is the suite-level refactor of the current Merlino/MORPHEUS codebase.
The first releases keep Merlino compatibility aliases while the scientific
packages are separated behind stable service, CLI and manifest contracts.

## Initial Scope

- `oracle-core`: workspace layout, configuration, manifests, checksums and
  common errors.
- `oracle-chem`: atoms, masses, geometry, topology, rings and symmetry.
- `oracle-gicforge`: non-redundant GIC construction, frozen schemas, SYCART and
  B-matrix evaluation.
- `oracle-morpheus`: semiexperimental geometry refinement, constraints,
  predicates, parameter classes, diagnostics and ensemble refinement.
- `oracle-gf`: frozen-GIC GF/PED analysis.
- `oracle-gaussian`: Gaussian input/output adapters.
- `oracle-engines`: Fortran and vendored backend discovery/build wrappers.
- `oracle-dvr`: scan/grid to DVR workflows.
- `oracle-vpt2-vci`: QFF, VPT2/VCI and Davidson workflows.
- `oracle-gui`: GUI controllers and views only.

## Migration Rule

Scientific behavior is migrated package by package. Existing `merlino_*`
imports stay valid until ORACLE-native tests cover the new public APIs.

## Workspace Contract

Project workspaces use:

```text
inputs/
runs/
outputs/
reports/
cache/
logs/
```

Every workflow run writes an `oracle.run.v1` manifest containing input/output
paths, SHA256 hashes, parameters, backend metadata and status.

