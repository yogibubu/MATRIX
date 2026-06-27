# MATRIX Refactor Plan

Date: 2026-06-26

## Goal

Create a clean MATRIX repository from the legacy workspace without losing the
working scientific code, manuals, regression fixtures or GUI workflows.

MATRIX (**Molecular Analysis Toolkit for Reusable Integrated eXperiments**) is
the framework/package family. ORACLE remains the GUI/user-facing application
(**Operator for Routing, Analysis, Control, Launch and Exploration**), NEO is
the GIC coordinate engine (**Nonredundant Equivariant Orthogonalizer**) and LINK
is the preprocessing/import layer (**Loader for Interoperable Normalized
Knowledge**). Existing module names remain compatibility aliases until the
package boundaries are stable and covered by tests. See
`ADR-0009-MATRIX-NAMING-TRANSITION.md`.

## Current Situation

The repository already contains useful MATRIX package boundaries:

- `matrix_core`: configuration, workspace layout, manifests and shared errors.
- `matrix_neo`: GICForge service, frozen GIC schemas and B-matrix evaluation.
- `matrix_morpheus`: MORPHEUS/SEfit single- and multi-structure refinement.
- `matrix_gf`: frozen-GIC GF/PED analysis.
- `matrix_qm`: shared QM tensor sections for Hessians, normal modes and QFF.
- `matrix_gaussian`: Gaussian adapters.
- `matrix_engines`: backend discovery and wrappers.
- `matrix_dvr`: DVR workflow wrapper.
- `matrix_vpt2_vci`: QFF, VPT2/VCI and Davidson prototypes.
- `matrix_oracle`: newer workflow dashboard.

The same tree also still contains legacy or mixed-responsibility areas:

- `gui` and `advanced` still orchestrate scientific services directly.
- `geometry`, `topology` and the legacy fitting stack overlap and create
  cross-package dependencies.
- `working`, `tmp`, generated paper artifacts and LaTeX build files pollute
  the checkout.
- `doc/papers` mixes source manuscripts, generated tables, analysis results and
  build products.
- `puckering_dvr` is a large vendored backend and should be treated as an
  engine, not as ordinary ORACLE application code.

The current worktree is dirty. Before moving files, preserve it as a historical
state and separate code changes from generated outputs.

## Manuals Read

The `newmsr_overleaf` manuals define the target scientific contracts:

- MORPHEUS Manual: single-structure semiexperimental refinement, constraints,
  predicates, parameter classes, GIC and symmetry-Cartesian coordinate models,
  diagnostics and CLI/GUI workflow.
- GICForge Manual: deterministic non-redundant GIC construction, ring and
  butterfly coordinates, symmetry adaptation, SYCART and the Python/Fortran
  identity contract.
- GF/PED Manual: Cartesian Hessian plus frozen GIC definition, Pulay-style
  scaling, Wilson GF/PED reports and CSV outputs.
- Multi-Structure MORPHEUS Manual: shared class-correction refinement across
  related molecules or conformers, priors, hard constraints, synthon `Zeff`
  typing, ensemble diagnostics and paper artifact generation.

These manuals should become ORACLE contract documentation, not side artifacts.

## Proposed ORACLE Repository Layout

Use a monorepo first. Split into multiple repositories only after interfaces,
tests and release boundaries are stable.

```text
ORACLE/
  pyproject.toml
  README.md
  scripts/
    matrix_env.sh
  docs/
    architecture/
    manuals/
    papers/
    archive/oracle3/
  packages/
    matrix-core/
      src/matrix_core/
    matrix-chem/
      src/matrix_chem/
    matrix-neo/
      src/matrix_neo/
    matrix-morpheus/
      src/matrix_morpheus/
    matrix-gf/
      src/matrix_gf/
    matrix-qm/
      src/matrix_qm/
    matrix-gaussian/
      src/matrix_gaussian/
    matrix-molpro/
      src/matrix_molpro/
    matrix-mrcc/
      src/matrix_mrcc/
    matrix-fragments/
      src/matrix_fragments/
    matrix-rovib/
      src/matrix_rovib/
    matrix-thermo/
      src/matrix_thermo/
    matrix-engines/
      src/matrix_engines/
    matrix-dvr/
      src/matrix_dvr/
    matrix-vpt2-vci/
      src/matrix_vpt2_vci/
    matrix-oracle/
      src/matrix_oracle/
  engines/
    fortran/gicforge/
    fortran/dvr/
    fortran/vpt2_vci/
    vendored/puckering_dvr/
  examples/
    morpheus/
    gicforge/
    gf/
    dvr/
    vpt2_vci/
  benchmarks/
    semiexp_msr/
  tests/
    fixtures/test_molecules/
    integration/
    regression/
  tools/
    oracle_run.py
    migrate_oracle.py
```

## Central Architecture Constraint

All MATRIX modules must reuse the same libraries for the same tasks. The suite
must not grow duplicate XYZ parsers, topology builders, isotope tables, GIC
builders, Gaussian parsers, manifest writers or backend launchers.

This does not forbid intentional dual scientific kernels. Some modules may have
both Python and strict Fortran77 implementations. In ORACLE these are treated as
backends behind one service contract, not as separate tools with separate data
models. They must consume the same enriched XYZ sections, use the same shared
libraries around the kernel, record backend metadata in the manifest and share
identity/regression tests.

The canonical communication file is an enriched XYZ container. It starts with a
plain XYZ block and is progressively enriched by named uppercase sections. Each
tool owns only its own section, replaces only that section and preserves all
others. External formats are import/export adapters; downstream MATRIX modules
consume the enriched XYZ container.

The implemented standalone tool contracts are recorded in
`matrix_core.tool_contracts`. That registry is intentionally small and
machine-readable: it lists current package names, future names where decided,
required/optional sections, produced sections and the canonical standalone CLI
entry point for each tool.

This constraint is formalized in
`ADR-0001-SHARED-LIBRARIES-AND-XYZ-CONTAINER.md` and
`MATRIX_XYZIN_CONTAINER.md`.

The preprocessing layer is LINK. It normalizes external sources into an
Avogadro-compatible enriched XYZ, determines symmetry with explicit thresholds,
builds topology once, saves synthons/descriptors and leaves downstream tools to
reuse the saved sections.

## Subproject Responsibilities

`matrix-core`

- Workspace layout: `inputs/`, `runs/`, `outputs/`, `reports/`, `cache`,
  `logs`.
- Run manifests, checksums, config, logging and shared typed errors.
- Shared sectioned enriched-XYZ utilities.
- No chemistry, GUI, Gaussian or Fortran-specific logic.

`matrix-chem`

- Atoms, masses, isotopes, XYZ, inertia, topology, rings, symmetry and primitive
  geometry utilities.
- Owns canonical geometry data models and plain/enriched XYZ parsers.
- Owns the unified Z-matrix parser. GUI, Gaussian and legacy adapters must call
  this parser rather than keeping private Z-matrix readers.
- This should absorb the stable parts of `geometry`, `topology` and topology
  pieces now under `matrix_morpheus`.
- Descriptor topology and atomic synthons are first-class chemistry primitives,
  not side tools. `AtomicSynthons`, descriptor parameters, aromaticity and
  fragment/synthon signatures belong here so MORPHEUS, GICForge, fragment
  search, GUI diagnostics and future ML/data workflows reuse the same
  descriptors.
- Owns `#VALIDATION`, the post-preprocessing gate that checks the enriched XYZ
  state before GICForge and later tools consume the molecule.
- `topology_reporting` is migrated as source material but should not become a
  public API until the remaining `survibfit` reporting dependencies are moved
  or replaced by MATRIX services.

`matrix-neo`

- GIC construction, ring numbering, symmetry labels, SYCART, frozen schemas,
  B-matrix evaluation and Python/Fortran comparison.
- Owns the GICForge public service, normalized schema and Python/Fortran77
  backend contract.
- Starts only from an enriched XYZ with `#VALIDATION STATUS PASS`; it then
  writes frozen `#GIC` and optional `#SYCART` sections for downstream tools.
- May create Gaussian input on request, but the file-format serialization is
  delegated to `matrix-gaussian`.

`matrix-morpheus`

- Single-structure and multi-structure semiexperimental refinement.
- Consumes frozen coordinate models from `matrix-neo` or
  symmetry-Cartesian data.
- Owns constraints, predicates, parameter classes, robust least squares,
  diagnostics, reports and reference-library search.

`matrix-gf`

- Frozen-GIC Hessian transformation, Wilson GF/PED, Pulay scaling, reports and
  CSV tables.
- Must stay independent from VPT2/VCI.

`matrix-qm`

- Shared enriched-XYZ contracts for Cartesian Hessians, normal modes and QFF
  data promoted by external QM adapters.
- Owns `#CARTESIAN_HESSIAN`, `#NORMAL_MODES` and `#QFF`.
- Scientific tools consume these sections instead of reparsing QM output files.

`matrix-gaussian`

- Gaussian input writing, log/FCHK/QFF parsing and normalized correction tables.
- Gaussian is an adapter, not an internal model.
- Gaussian geometry/log parsers must return shared MATRIX models such as
  `matrix_chem.MolecularGeometry`; downstream tools must not parse Gaussian text
  directly.
- Gaussian Z-matrix input delegates to the shared `matrix-chem` Z-matrix parser.
- Gaussian input generation from GICForge consumes enriched XYZ plus `#GIC`.
- Gaussian rovibrational log data is promoted to shared `#VIBRATIONAL`,
  `#ROTATIONAL` and `#DELTABVIB` sections before GF, Thermo, SEfit/MORPHEUS or
  anharmonic workflows consume it.

`matrix-molpro` / `matrix-mrcc`

- Program-specific output adapters for Molpro and MRCC geometry import.
- Return shared `matrix_chem.MolecularGeometry` objects and feed LINK
  preprocessing, which writes `#BASIC`, `#SOURCE`, `#SYMMETRY`, `#TOPOLOGY` and
  `#SYNTHONS`.
- Downstream GF, SEfit/MORPHEUS and anharmonic workflows must consume enriched
  XYZ sections, not Molpro/MRCC text.

`matrix-fragments`

- Nano-LEGO-like fragmentation, fragment-library search and future assembly
  workflows.
- It is a client of `#TOPOLOGY` and `#SYNTHONS`, not a separate topology
  engine.
- LCB25 molecules are first imported/preprocessed by LINK, then indexed
  here as full references or topology/synthon fragments.
- Query molecule fragmentation and LCB25 fragment-library lookup are inverse
  uses of the same service contract.

`matrix-rovib`

- Rotational and vibrational section contracts, DeltaVib/alpha bridge values,
  Coriolis and Q-cent compatibility.
- Provides standalone `xyzin` readers/writers for `#ROTATIONAL` and
  `#VIBRATIONAL`, plus summary tooling for existing MATRIX containers.
- Migrates ORACLE `geometry/rotational_pipeline.py`, `vibrational.py`,
  `vib_anh.py`, `rovib_pipeline.py`, `coriolis.py` and `qcent.py`.
- Treats external CeDiTT/alpha-resonances payloads as imported data that enrich
  `#ROTATIONAL`, not as private parser logic in downstream tools.

`matrix-thermo`

- Thermochemistry from enriched XYZ state.
- Reads `#BASIC`, `#ROTATIONAL` and optional `#VIBRATIONAL`.
- Owns `#THERMO` and migrates ORACLE `geometry/thermo_*` modules.

`matrix-engines`

- Discovery, build checks and subprocess wrappers for Fortran and vendored
  engines.
- GUI and scientific packages should call this layer instead of executing
  binaries directly.
- Backend metadata and executable/source checksums are recorded in manifests.

`matrix-dvr`

- Scan/grid to DVR workflow, Cremer-Pople labeling and output readers.
- Treat `puckering_dvr` as a vendored engine or adapter backend.

`matrix-vpt2-vci`

- Normal-mode QFF, VPT2/VCI basis management and Davidson diagonalization.
- Consumes normalized force-field data from `matrix-gaussian`.

`matrix-oracle`

- PySide6 controllers and views only.
- No scientific algorithms, no direct Fortran calls, no duplicate parsers.

## Pipeline Contract

Every project workflow should follow the same shape:

```text
inputs/
  user-provided geometries, job files, Gaussian logs, observation tables
runs/
  timestamped immutable run directories
outputs/
  canonical machine-readable outputs
reports/
  human-readable text, HTML, LaTeX or PDF summaries
cache/
  reusable derived data
logs/
  application and backend logs
```

Every run writes an `oracle.run.v1` manifest with:

- workflow name and schema version
- input paths and SHA256 hashes
- output paths and SHA256 hashes
- parameters
- backend executable/source metadata
- git commit and dirty flag when available
- status and user-facing messages

During compatibility migration, accept and emit `oracle.run.v1` where needed,
but the new ORACLE manifest should be the forward contract.

Every workflow should also accept or produce an enriched XYZ file as the
canonical state handoff. Workflow-specific files may be generated in `runs/`,
`outputs/` or `reports/`, but the reusable state must be written back to the
appropriate XYZ section through the shared section API.

Standalone mode is mandatory for scientific packages. SEFit/MORPHEUS, GF/PED,
Thermo, DVR and VPT2/VCI must accept a sufficiently populated `xyzin` directly
and validate only the sections they consume; a full MATRIX project workspace is
recommended but not required for those direct runs.

## Migration Phases

### Phase 0: Freeze and Inventory

1. Create a clean archive or branch from the current dirty `oracle-legacy` state.
2. Split pending changes into groups: code, examples, benchmarks, papers,
   generated artifacts and local runtime outputs.
3. Add missing ignore rules for LaTeX build files and runtime outputs that are
   currently showing up in `git status`.
4. Copy the manuals from `newmsr_overleaf` into `doc/manuals` or confirm that
   the existing `doc/manuals` copy is canonical.

Exit criteria:

- Historical state preserved.
- No generated runtime output mixed with source changes.
- Existing `./freeze_check.sh` behavior documented, even if not yet green.

### Phase 1: ORACLE Skeleton

1. Create the ORACLE monorepo skeleton with workspace-aware `pyproject.toml`.
2. Add package placeholders and compatibility aliases from legacy imports to
   MATRIX packages.
3. Move only docs and metadata first; do not move scientific code in the first
   commit.
4. Add `python -m matrix` as the initial ORACLE CLI.

Exit criteria:

- `python -m matrix --help` works.
- Existing compatibility imports still work.
- No scientific behavior changed.

### Phase 2: Core and Engines

1. Extract core services into `matrix-core`.
2. Extract engine discovery and backend build checks into `matrix-engines`.
3. Normalize manifest schema to `oracle.run.v1` while preserving legacy
   manifest readers.
4. Replace direct subprocess calls in services with engine wrappers.

Exit criteria:

- Core tests pass.
- GICForge/DVR/Fortran build discovery tests pass.
- Manifest compatibility test covers both schemas.

### Phase 3: Chemistry Foundation

1. Merge stable `geometry`, `topology` and selected `matrix_morpheus.topology`
   services into `matrix-chem`.
2. Remove circular dependencies between topology and `matrix_morpheus`.
3. Define one public molecular model for atoms, coordinates, masses, graph,
   rings and symmetry metadata.

Exit criteria:

- `matrix-neo` and `matrix-morpheus` can consume `matrix-chem`.
- No package imports legacy `geometry` or `topology` directly except
  compatibility wrappers.

### Phase 4: GICForge and MORPHEUS

1. Move GICForge services and frozen schema models into `matrix-neo`.
2. Move MORPHEUS single-structure and ensemble workflows into
   `matrix-morpheus`.
3. Keep the manuals as contract docs and add examples that match each manual.
4. Add regression fixtures for GICForge, SYCART, single MORPHEUS and ensemble
   MORPHEUS. The imported ORACLE `test_molecules` corpus is the first
   demanding GICForge/parser fixture set.

Exit criteria:

- `matrix gicforge define`, `matrix morpheus fit` and
  `matrix morpheus ensemble` run without importing GUI modules.
- Existing semiexp tests pass under ORACLE names and ORACLE aliases.

### Phase 5: Analysis Engines

1. Move GF/PED into `matrix-gf`.
2. Move Gaussian adapters into `matrix-gaussian`.
3. Move DVR orchestration into `matrix-dvr`.
4. Move VPT2/VCI into `matrix-vpt2-vci`.
5. Keep GF independent from VPT2/VCI and keep Gaussian as an adapter.

Exit criteria:

- `matrix gf`, `matrix gaussian summary`, `matrix dvr`, `matrix vpt2-vci`
  expose non-GUI CLI workflows.
- Each writes an ORACLE manifest and has small fixture tests.

### Phase 6: GUI Rewire

1. Replace legacy `gui` and `advanced` direct logic with calls into
   `oracle-*` services.
2. Keep GUI state, view models and file selection in `matrix-oracle`.
3. Remove scientific parsing, fitting and direct Fortran execution from GUI
   classes.

Exit criteria:

- GUI smoke tests pass.
- Service/CLI tests cover the scientific behavior behind each GUI action.

### Phase 7: Cleanup and Release

1. Move old manuscripts and historical files into `docs/archive/oracle3`.
2. Remove duplicated generated artifacts from source control.
3. Keep benchmark inputs and golden outputs, but regenerate reports through
   explicit commands.
4. Publish ORACLE 0.1 as a compatibility release with ORACLE aliases.

Exit criteria:

- Clean `git status` after a documented build/test run.
- New clone can run tests without local `working/` state.
- Documentation explains legacy-to-ORACLE mapping.

## Immediate First Tasks

1. Commit or archive the current dirty state before any moves.
2. Add a repository hygiene patch:
   ignore LaTeX build products, keep generated paper outputs out of ordinary
   status, and untrack runtime-only files after review.
3. Create ORACLE skeleton in a separate directory or branch.
4. Copy transition plans, package architecture notes, repository layout notes
   and manuals into ORACLE docs.
5. Implement only `matrix-core` and compatibility aliases first.
6. Run the smallest validation set before each migration step.

## Risk Register

- Dirty worktree: high risk of losing active scientific edits if cleanup and
  migration are mixed.
- GUI coupling: medium-high risk; GUI imports many scientific packages and must
  be rewired last.
- Topology/GIC coupling: high risk; `matrix_neo` and `matrix_morpheus` still
  depend on `matrix_morpheus` and legacy topology code.
- Generated paper artifacts: medium risk; useful for publications but noisy for
  source control.
- Fortran/vendored engines: medium risk; treat as stable engines with narrow
  wrappers.
- Rename risk: high if modules are renamed too early. Keep ORACLE aliases until
  ORACLE tests are green.

## Validation Matrix

- Core: workspace layout, config, manifest checksum tests.
- Chemistry: XYZ, masses, topology, ring numbering and symmetry tests.
- GICForge: Python/Fortran identity, SYCART, frozen schema and B-matrix tests.
- MORPHEUS: single-structure fit, constraints, predicates, class constraints,
  robust loss, leave-one-out and ensemble tests.
- GF/PED: frozen-GIC FCHK workflow, scaling and CSV export tests.
- Gaussian: input writer, log/FCHK/QFF parser tests.
- DVR: Gaussian log/grid to DVR args, manifest and output reader tests.
- VPT2/VCI: QFF load, basis construction, VPT2, VCI and Davidson tests.
- GUI: smoke tests only after service/CLI tests pass.
