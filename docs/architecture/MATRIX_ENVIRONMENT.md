# MATRIX Environment Helpers

MATRIX keeps sourceable shell helpers in the repository instead of editing a
personal shell startup file implicitly. Source the helper once per shell:

```bash
source /Users/vincenzobarone/MATRIX/scripts/matrix_env.sh
matrix-set
```

The helper defines:

- `matrix-set`: create `MATRIX_VENV` when missing, activate the first available
  MATRIX environment, install runtime dependencies including RDKit when needed,
  export `MATRIX_HOME`, add package `src` directories to `PYTHONPATH`, add local
  tools to `PATH`, and `cd` to the repo.
- `matrix-unset`: deactivate the MATRIX virtualenv/conda env when active,
  restore the previous `PATH`, `PYTHONPATH` and working directory saved by
  `matrix-set`.
- `matrix-run`: launch `matrix_oracle.app` when it exists, otherwise dispatch to
  `python -m matrix`.
- `matrix-cli`: run `python -m matrix` directly.
- `matrix-run-bg`: launch the same target in background with a log file.
- `matrix-run-check`: verify the current Python can import the core scientific
  runtime stack and report optional external viewers such as Molden.
- `matrix-install-runtime-deps`: explicitly install/upgrade the core runtime
  stack in the active MATRIX environment.
- `matrix-test`, `matrix-test-all`: run focused or full tests in the activated
  environment.
- `matrix-clean`: remove Python cache by default; with `outputs` or `all`, also
  remove ordinary runtime logs while preserving `tests/fixtures`.
- `matrix-gic-corpus-list`: list imported demanding GIC input files.

Default variables use MATRIX names; the script still maps old `ORACLE_*`
variables to the new names during the transition:

```bash
MATRIX_HOME=/Users/vincenzobarone/MATRIX
MATRIX_VENV=$HOME/.venvs/matrix
MATRIX_CONDA_ENV=matrix_26
MATRIX_PYTHON=python3
MATRIX_AUTO_CREATE_VENV=1
MATRIX_AUTO_INSTALL_RUNTIME_DEPS=1
MATRIX_AUTO_INSTALL_GUI_DEPS=0
MATRIX_RUNTIME_DEPS="numpy scipy matplotlib pandas sympy pytest rdkit"
MATRIX_GUI_DEPS="PySide6 pytest-qt"
```

Dense Hermitian diagonalization is routed through
`matrix_core.diagonalizer`. The default policy uses SciPy/NumPy on CPU and can
use CuPy or PyTorch GPU backends when they are already installed in the active
environment. GPU packages are not installed by `matrix-set` because the correct
wheel depends on the local hardware and driver stack.

Useful controls:

```bash
MATRIX_DIAGONALIZER_BACKEND=auto
MATRIX_DIAGONALIZER_GPU_MIN_SIZE=128
MATRIX_DIAGONALIZER_STRICT_GPU=0
```

Molden is an optional external viewer for the Electronic Spectroscopy GUI. On
macOS it needs a Molden executable in `PATH` and XQuartz for the X11 display.
`matrix-run-check` reports both conditions, but it does not install XQuartz
because the macOS package installer requires an interactive administrator
password.

Set `MATRIX_AUTO_CREATE_VENV=0` to prevent `matrix-set` from creating a
virtualenv. Set `MATRIX_AUTO_INSTALL_RUNTIME_DEPS=0` to prevent automatic
runtime dependency installation. RDKit is part of the runtime dependency set
because LINK uses it for SMILES imports.

Set `MATRIX_AUTO_INSTALL_GUI_DEPS=1` before `matrix-set` only when automatic
GUI dependency installation is desired. For routine development, explicit
dependency installation is safer:

```bash
matrix-install-runtime-deps
matrix-install-gui-deps
```

To make the commands permanent in bash, add only this line to `~/.bashrc`:

```bash
source /Users/vincenzobarone/MATRIX/scripts/matrix_env.sh
```

For zsh, put the same line in `~/.zshrc`.

Legacy `oracle-set`, `oracle-run`, `oracle-test-all` and `scripts/oracle_env.sh`
remain wrappers for existing shells, but new documentation should use the
MATRIX names.
