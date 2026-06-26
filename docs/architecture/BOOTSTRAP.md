# Bootstrap

This skeleton is intentionally light. It should prove package boundaries before
scientific code is moved.

## Local Smoke Test

```bash
cd /Users/vincenzobarone/ORACLE
PYTHONPATH=packages/oracle-core/src python -m pytest
```

## Create A Workspace

```bash
cd /Users/vincenzobarone/ORACLE
PYTHONPATH=packages/oracle-core/src python tools/oracle_run.py init /tmp/oracle-demo
```

## Delegate To Merlino During Migration

When the Merlino repository is importable in `PYTHONPATH`, the temporary CLI can
delegate:

```bash
PYTHONPATH=/Users/vincenzobarone/merlino3.0 python tools/oracle_run.py merlino --help
```

