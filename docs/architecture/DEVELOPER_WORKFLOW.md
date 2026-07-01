# MATRIX Developer Workflow

Every new workflow should be added in this order:

1. Define a non-Qt service API with dataclasses for normalized inputs/outputs.
2. Add CLI access under `python -m matrix`.
3. Write an `oracle.run.v1` manifest with input/output checksums.
4. Add focused tests using small fixtures.
5. Add GUI wiring only after the service and CLI are tested.
6. Document file contracts and benchmark commands.
7. Run `matrix-test-all` or `python -m pytest` before committing.

Use the sourceable shell helpers when working interactively:

```bash
source /Users/vincenzobarone/MATRIX/scripts/matrix_env.sh
matrix-set
matrix-run --help
matrix-test-all
```

Rules:

- GUI classes must not contain scientific algorithms.
- Gaussian, Fortran and other external formats are adapters, not internal data
  models.
- Coordinate work must respect `docs/architecture/INTERNAL_COORDINATE_FRAMEWORK.md`:
  topology, geometry, coordinate generation, reduction, symmetrization and
  downstream analysis are separate stages.
- Prefer typed ORACLE errors from `matrix_core.errors`.
- Important numerical fits should use shared MATRIX numerical primitives for
  damped normal equations, step limiting and rank/condition diagnostics unless
  a workflow has a documented reason to use a specialized solver.
- Run manifests must be built through `matrix_core.build_run_manifest`. The
  shared builder records MATRIX version, command line, Python version,
  input/output checksums, detected `xyzin` sections and external backend
  metadata.
- Use `tests/fixtures/test_molecules` for demanding GICForge/parser
  regressions instead of inventing ad hoc molecule inputs.
- Store new project outputs under `inputs/`, `runs/`, `outputs/`, `reports/`,
  `cache/` or `logs/`.

## Internal-Coordinate Change Checklist

Before changing NEO/GICForge, GF coordinate handling, topology clients or
symmetry code, answer these questions in the commit, issue or review notes.

- Does the change preserve Merlino/GICForge numerical behaviour unless the
  requested scientific change explicitly requires otherwise?
- Is there already a module that owns this responsibility?
- Are topology, geometry, coordinate generation, numerical reduction and
  symmetrization still separated?
- Is the implementation reusing existing parsers/generators/policies rather
  than duplicating them?
- If a coordinate type changed, is there one clear logical generator and one
  registry entry for it?
- If behaviour changed, is there a focused regression test and, when relevant,
  a Python/Fortran comparison?
- Does the frozen xyzin contract still contain enough metadata for downstream
  tools to work without reparsing external program input?
- For symmetry changes, do the tests verify labels, total-symmetric selection
  and absence of forbidden cross-irrep GF couplings where applicable?
- For large-amplitude/local-mode changes, do the tests verify both F and G
  coupling diagnostics?

## Continuous Integration

The GitHub Actions workflow `.github/workflows/ci.yml` is the reference
repository health check. It performs the same bootstrap expected from a fresh
checkout:

```bash
source scripts/matrix_env.sh
matrix-set
python -m matrix --help
python -m oracle --help
matrix-run --help
python -m ruff check packages tests
python -m pytest -q
```

The CI job enables `MATRIX_AUTO_INSTALL_GUI_DEPS=1`, `QT_QPA_PLATFORM=offscreen`
and `MPLBACKEND=Agg`, installs `gfortran` and Qt headless system libraries, and
therefore covers CLI imports, ORACLE GUI importability, the ruff correctness
gate and the Fortran-backed regression tests when the vendored sources can be
compiled.
