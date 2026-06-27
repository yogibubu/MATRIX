# Bootstrap

This skeleton is intentionally light. It should prove package boundaries before
scientific code is moved.

## Local Smoke Test

```bash
cd /Users/vincenzobarone/MATRIX
PYTHONPATH=packages/matrix-core/src python -m pytest
```

## Shell Helpers

ORACLE provides sourceable shell helpers for environment setup, launch checks
and tests:

```bash
source /Users/vincenzobarone/MATRIX/scripts/matrix_env.sh
matrix-set
matrix-run --help
matrix-run-check
matrix-test-all
```

See `docs/architecture/MATRIX_ENVIRONMENT.md`.

## Create A Workspace

```bash
cd /Users/vincenzobarone/MATRIX
PYTHONPATH=packages/matrix-core/src python -m matrix init /tmp/oracle-demo
```

## Inspect The CLI

The temporary CLI is available directly from the repository:

```bash
python -m matrix --help
```
