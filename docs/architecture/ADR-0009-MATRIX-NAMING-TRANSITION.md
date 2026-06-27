# ADR-0009: MATRIX Framework Naming Transition

Date: 2026-06-27

## Status

Accepted and implemented. Physical MATRIX package renames are active, with
ORACLE-era entry points kept as compatibility aliases.

## Context

The repository originally used ORACLE-oriented names for the framework,
runtime packages and GUI. The refactory now separates the framework,
user-facing GUI and scientific tools while preserving compatibility for
existing scripts and xyzin workflows.

## Decision

The public naming model is:

- **MATRIX** becomes the framework/package family name:
  **Molecular Analysis Toolkit for Reusable Integrated eXperiments**.
- **ORACLE** remains the user-facing GUI/application name:
  **Operator for Routing, Analysis, Control, Launch and Exploration**.
- **GICForge** is renamed to **NEO**:
  **Nonredundant Equivariant Orthogonalizer**.
- **ORACLE-Babel** is renamed to **LINK**:
  **Loader for Interoperable Normalized Knowledge**.
- Other tools should use Matrix-saga character names when a rename is useful
  and does not obscure the scientific contract.

NEO is assigned to GICForge rather than to the GUI because that tool is the
coordinate engine that builds, projects, symmetrizes and reduces the molecular
internal-coordinate representation. ORACLE is a better GUI name because the GUI
routes users through project state, launches tools and exposes the scientific
state without owning parser or kernel logic.

LINK replaces "Babel" because "Babel" describes translation but does not belong
to the Matrix naming universe. LINK is the MATRIX entry adapter: it connects
external XYZ, QM outputs, Z-matrices, SMILES/RDKit and databases such as LCB25
to the normalized xyzin molecular state.

Runtime packages have been physically renamed from `oracle-*` to `matrix-*`.
Python imports now use `matrix_*` modules as the primary API. Legacy
`oracle_*` modules remain as thin compatibility shims and should not receive
new implementation code.

The active command layer is:

- `matrix` is the primary framework CLI.
- `matrix neo ...` is an alias for `matrix gicforge ...`.
- `neo ...` is an installable console alias for the GICForge/NEO coordinate
  tool.
- `matrix link ...` is the documented primary command for LINK.
- `matrix babel ...` and `oracle babel ...` remain compatibility commands for
  old workflows.
- `matrix-set`, `matrix-run`, `matrix-cli` and `matrix-test-all` are the
  primary shell helpers; `oracle-set`, `oracle-run`, `oracle-cli` and
  `oracle-test-all` remain wrappers.

`matrix_core.tool_contracts` records current names, planned names and
compatibility aliases after the physical package renames. It is the source of
truth for tool ownership and migration planning.

## Consequences

- Public contracts remain stable during scientific porting.
- Renames happen as a planned compatibility migration, not as incidental churn.
- Docs can introduce tools as, for example, "NEO, formerly GICForge", with
  compatibility aliases kept during the transition.
- Docs can introduce "LINK, formerly ORACLE-Babel"; the old Babel command
  names stay available for scripted workflows.
- The GUI may be branded ORACLE while running on the MATRIX framework.
- New code should import `matrix_*` modules and use MATRIX command names.
- Compatibility shims prevent immediate breakage of downstream ORACLE scripts.
