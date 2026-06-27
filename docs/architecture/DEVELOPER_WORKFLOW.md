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
- Prefer typed ORACLE errors from `matrix_core.errors`.
- Important numerical fits should use shared MATRIX numerical primitives for
  damped normal equations, step limiting and rank/condition diagnostics unless
  a workflow has a documented reason to use a specialized solver.
- Use `tests/fixtures/test_molecules` for demanding GICForge/parser
  regressions instead of inventing ad hoc molecule inputs.
- Store new project outputs under `inputs/`, `runs/`, `outputs/`, `reports/`,
  `cache/` or `logs/`.

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
python -m pytest -q
```

The CI job enables `MATRIX_AUTO_INSTALL_GUI_DEPS=1`, `QT_QPA_PLATFORM=offscreen`
and `MPLBACKEND=Agg`, installs `gfortran` and Qt headless system libraries, and
therefore covers CLI imports, ORACLE GUI importability and the Fortran-backed
regression tests when the vendored sources can be compiled.
