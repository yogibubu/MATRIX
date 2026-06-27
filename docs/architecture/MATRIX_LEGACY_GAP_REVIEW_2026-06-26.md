# MATRIX Legacy Gap Review

Date: 2026-06-26

Source reviewed: `/Users/vincenzobarone/oracle-legacy`

## Summary

The main package families are represented in ORACLE, but the current ORACLE
plan under-described several ORACLE responsibilities:

- standalone scientific runs from an existing `xyzin`;
- rovibrational/DeltaVib compatibility utilities;
- thermochemistry;
- local semiexperimental reference libraries;
- paper/benchmark artifact generation;
- Gaussian job launching/status helpers.

These are now tracked as explicit architecture responsibilities. The Gaussian,
Thermo and Rovib items have MATRIX-native ports; remaining items should migrate
after the shared parser, topology and section contracts are stable.

## Reviewed ORACLE Areas

Core/package areas:

- `matrix_core`: CLI, config, manifests, workspace, `xyzin` sections and
  geometry helpers.
- `matrix_neo`: GICForge Python/Fortran service, frozen GIC schema, symmetry,
  B-matrix evaluation.
- `matrix_morpheus`: SEFit/MORPHEUS, ensemble fitting, MSR import, Kraitchman,
  reference library search, paper benchmarks.
- `matrix_gf`: HessianInput, Wilson GF/PED, internal-coordinate transform,
  reports and CSV export.
- `matrix_vpt2_vci`: normalized QFF consumers, VPT2, VCI, Davidson and workflow.
- `matrix_dvr`: DVR request and manifest preparation.
- `matrix_gaussian`: Gaussian input/log/FCHK/QFF adapters plus job
  discovery/status/run helpers.
- `matrix_engines`: backend discovery and build checks.
- `geometry`: rotational, vibrational, rovib, Coriolis, Q-cent and
  thermochemistry utilities.
- `topology` and `matrix_morpheus/topology`: graph/rings/synthons/descriptors.
- `advanced` and `gui`: GUI orchestration, launchers and viewers.

Data/docs:

- `data/se_geometries`: local SE reference geometry library.
- `benchmarks/semiexp_msr`: regression and paper benchmark inputs/goldens.
- `doc/XYZIN_FORMAT.md`: historical `xyzin` section contract.
- `doc/ROVIB_COMPATIBILITY_INTERFACE_2026-03-24.md`: DeltaVib/alpha bridge
  contract.
- `doc/VPT2_VCI_CORE.md`: GF/VPT2/VCI separation and adapter policy.
- `doc/FORTRAN_BACKENDS.md`: active Fortran backends and Python ownership.

## Gaps Added To ORACLE Plan

`matrix-rovib`

- Owns normalized `#ROTATIONAL`, `#VIBRATIONAL`, `#DELTABVIB`, `#CORIOLIS`
  and `#QCENT` contracts.
- Migrates ORACLE `geometry/rotational_pipeline.py`, `vibrational.py`,
  `vib_anh.py`, `rovib_pipeline.py`, `coriolis.py` and `qcent.py`.
- Current ORACLE port includes Merlino-compatible `vibin` generation from the
  shared Gaussian FCHK adapter, sparse Coriolis terms, QCENT Watson constants,
  direct vibrational DOS and rovibrational DOS/Q(T) commands.
- Treats CeDiTT/alpha-resonances output as an external compatibility payload,
  not as ORACLE-owned theory unless explicitly migrated later.

`matrix-thermo`

- Owns `#THERMO`.
- Migrates `geometry/thermo_trasl.py`, `thermo_rot.py`, `thermo_vib.py`,
  `thermo_pipeline.py` and `thermo_writer.py`.
- Reads `#BASIC`, `#ROTATIONAL` and optional `#VIBRATIONAL` from `xyzin`.
- Current ORACLE port keeps Merlino's translational, rotational
  quantum/classical crossover and harmonic vibrational formulas while using
  `matrix-core`, `matrix-chem` and `matrix-rovib` section contracts.

Standalone `xyzin` mode

- SEFit/MORPHEUS, GF/PED, Thermo, DVR and VPT2/VCI must be runnable from a
  sufficiently populated `xyzin` file.
- ORACLE pipeline commands create good state, but package CLIs must not require
  the user to rerun import/preprocess when the file already contains the needed
  sections.

Reference libraries

- LCB25 is already cached in ORACLE.
- ORACLE `data/se_geometries` remains a separate source to import or merge
  into the ORACLE library/index layer after topology/synthon validation.

Gaussian job helpers

- ORACLE has parsers, an input writer, job discovery/status/run helpers and
  `formchk` support in `matrix-gaussian`.
- FCHK/QFF parsing is owned by `matrix-gaussian`; `matrix_vpt2_vci.gaussian_qff`
  is retained only as a compatibility wrapper.

Paper/benchmark generation

- `semiexp-benchmark`, `semiexp-ensemble-paper` and generated publication
  artifacts are not core runtime services.
- They should migrate as benchmark/report commands after MORPHEUS core APIs are
  stable.

## Not Missing, But Still Pending

- Python/Fortran77 dual backends are already allowed by ADR-0002.
- GICForge, GF/PED, DVR and VPT2/VCI package placeholders already exist.
- Fragment/Nano-LEGO direction is covered by `matrix-fragments`.
- Unified geometry/Gaussian/Z-matrix parsing now includes Molpro and MRCC
  geometry/charge/multiplicity adapters. Other non-Gaussian Hessian,
  anharmonic and property adapters still need migration as source formats are
  identified.
- GIC symmetry is an ORACLE strength and must preserve Merlino3.0 behavior.
  The current MATRIX projector covers tested low-symmetry, `C/D/Dnh/Dnd`
  families plus matrix-classified `Td`, `Oh` and `Ih`. Remaining parity work is
  strict Fortran77 projector diagnostics equivalent to Python and regression
  fixtures for ring/butterfly/special coordinates under symmetry.

## Legacy Namespace Cleanup 2026-06-26

Closed in the first MORPHEUS cleanup pass:

- `merlino_core` service imports in MORPHEUS were replaced by `matrix_core`
  (`manifest`, `numerics`, `xyzin_geometry`, `isotopologues`, `paths`,
  scientific errors and BDPCS3 parameters).
- Generic geometry utilities used by SEfit/Kraitchman were moved to
  `matrix_chem` (`average_atomic_masses`, `isotopes_table`,
  `physical_constants`, `structure`, `inertia`, `rotational`).
- XYZ tuple I/O and XYZ writing now come from `matrix_chem.geometry_io` instead
  of `survibfit.modify_geom`.
- Topology/radius/element imports in MORPHEUS and the GIC kernel now use
  `matrix_chem.topology`.
- GICForge executable discovery now uses `matrix_engines` instead of
  `merlino_fortran`.
- Runtime aliases for `merlino_core`, `merlino_fortran`, `geometry` and
  `topology` were removed.

Validation for this cleanup pass: `PYTHONPATH=. pytest -q` returned
`125 passed`.

## Survibfit/GIC Runtime Port 2026-06-27

Closed in the GIC kernel cleanup pass:

- `survibfit` primitive, B-matrix, transform, local/global symmetry,
  similarity and supporting utilities were ported under
  `matrix_neo.survibfit`.
- The legacy-compatible GIC runtime service/model/symmetry postprocessor were
  ported under `matrix_neo.runtime`.
- MORPHEUS imports now target `matrix_neo.survibfit` and
  `matrix_neo.runtime` directly.
- The MORPHEUS runtime alias installer was removed; no `merlino_fit` or
  `merlino_gic` alias is installed at import time.
- The old MORPHEUS-vendored `legacy/merlino_fit/survibfit` and
  `legacy/merlino_gic` source trees were removed after the ORACLE port.
- The old MORPHEUS-vendored topology wrappers were removed; the active
  topology source is `matrix_chem.topology`.
- `matrix_neo.survibfit.pipeline` now uses `matrix_chem.topology`
  directly instead of loading a sibling topology package.
- Runtime GIC definitions now write `oracle.gic.definition.v1` while still
  accepting `merlino.gic.definition.v1` in readers for historical files.
- `MERLINO_FIT_*` runtime environment names were replaced by
  `ORACLE_GICFORGE_*`.

Remaining intentional compatibility:

- `LEGACY_*` schema constants with `merlino.*` values remain only where readers
  must accept historical job/checkpoint/GIC files.

Validation for this cleanup pass:

- `PYTHONPATH=. pytest -q tests/test_matrix_morpheus.py tests/test_matrix_neo.py tests/test_gic_regression_corpus.py`
  returned `46 passed`.
- `PYTHONPATH=. pytest -q` returned `125 passed`.

## VibInt and VarAnh Local Review 2026-06-27

Reviewed local sources:

- `/Users/vincenzobarone/vibint`
- `/Users/vincenzobarone/varanh`

`vibint` is a compact GF scaling utility centered on
`/Users/vincenzobarone/vibint/scale-gf.py`. It reconstructs a Gaussian/GDV
frequency-log GF problem from printed GIC, Wilson B, B-derivative, Hessian and
gradient blocks, then applies force-constant scaling and recomputes harmonic
frequencies. Its scientific core should not be ported into MATRIX because the
MATRIX GF path is already the stronger architecture:

- external Gaussian parsing is owned by the shared QM adapters, not by GF;
- NEO owns the frozen GIC definition and analytic B-matrix evaluation;
- GF consumes normalized `#GIC` and `#CARTESIAN_HESSIAN` sections;
- Pulay-style scaling already acts on GIC force constants and applies
  off-diagonal geometric-mean factors;
- `--local`, force thresholds and electrostatic/UFF-vdW subtraction are already
  implemented in MATRIX GF.

The useful `vibint` idea is user-facing rather than algorithmic: its rule
grammar for scaling by coordinate family, chemical subtype, wildcard and
priority is convenient. This should be considered as a MATRIX GF frontend that
translates typed rules into the existing `--scale`/`--scale-file` factors using
frozen `#GIC` labels, GIC names and `#SYNTHONS`. The Gaussian-log parser,
Tkinter GUI and private GF reconstruction in `vibint` should not be recycled.

`varanh` contains historical Fortran77 one-dimensional DVR and variational
anharmonic prototypes (`dvrHam.f`, `VCI.f`). The active MATRIX implementation
already covers the reusable parts in a cleaner way:

- `engines/fortran/dvr/path_dvr.f` provides the active downstream Fortran77 DVR
  kernel for nonperiodic sinc, periodic Fourier, distributed-Gaussian and 2D
  product-grid Hamiltonians;
- `packages/matrix-dvr` owns normalized `#DVR` state, manifests, direct runs
  and output collection;
- `engines/puckering_dvr/scripts/mw_path_dvr.py` handles Gaussian-derived
  mass-weighted paths, periodic/nonperiodic choices, 1D/2D grids and plotting;
- `packages/matrix-vpt2-vci` and `engines/fortran/vpt2_vci` already provide
  normalized VPT2/VCI basis generation, pruning, symmetry blocks and Davidson
  support.

No `varanh` code should be ported at this stage. It can remain a historical
source for formula checks or future golden regression cases, but it does not
identify a missing MATRIX DVR/VCI capability.
