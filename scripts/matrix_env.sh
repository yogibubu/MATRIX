# Source this file from bash or zsh to install MATRIX shell helpers.
#
# Example:
#   source /Users/vincenzobarone/MATRIX/scripts/matrix_env.sh
#   matrix-set

if [ -n "${BASH_SOURCE:-}" ]; then
    _matrix_env_file="${BASH_SOURCE[0]}"
elif [ -n "${ZSH_VERSION:-}" ]; then
    _matrix_env_file="${(%):-%x}"
else
    _matrix_env_file="$0"
fi

_matrix_env_dir="$(cd "$(dirname "$_matrix_env_file")/.." 2>/dev/null && pwd)"

export MATRIX_HOME="${MATRIX_HOME:-${ORACLE_HOME:-${_matrix_env_dir:-$HOME/MATRIX}}}"
export MATRIX_VENV="${MATRIX_VENV:-${ORACLE_VENV:-$HOME/.venvs/matrix}}"
export MATRIX_CONDA_ENV="${MATRIX_CONDA_ENV:-${ORACLE_CONDA_ENV:-matrix_26}}"
export MATRIX_PYTHON="${MATRIX_PYTHON:-${ORACLE_PYTHON:-python3}}"
export MATRIX_AUTO_CREATE_VENV="${MATRIX_AUTO_CREATE_VENV:-${ORACLE_AUTO_CREATE_VENV:-1}}"
export MATRIX_AUTO_INSTALL_RUNTIME_DEPS="${MATRIX_AUTO_INSTALL_RUNTIME_DEPS:-${ORACLE_AUTO_INSTALL_RUNTIME_DEPS:-1}}"
export MATRIX_AUTO_INSTALL_GUI_DEPS="${MATRIX_AUTO_INSTALL_GUI_DEPS:-${ORACLE_AUTO_INSTALL_GUI_DEPS:-0}}"
export MATRIX_RUNTIME_DEPS="${MATRIX_RUNTIME_DEPS:-${ORACLE_RUNTIME_DEPS:-numpy scipy matplotlib pandas sympy pytest rdkit}}"
export MATRIX_GUI_DEPS="${MATRIX_GUI_DEPS:-${ORACLE_GUI_DEPS:-PySide6 pytest-qt}}"

export ORACLE_HOME="${ORACLE_HOME:-$MATRIX_HOME}"
export ORACLE_VENV="${ORACLE_VENV:-$MATRIX_VENV}"
export ORACLE_CONDA_ENV="${ORACLE_CONDA_ENV:-$MATRIX_CONDA_ENV}"

matrix-package-path() {
    find "$MATRIX_HOME/packages" -mindepth 2 -maxdepth 2 -type d -name src 2>/dev/null \
        | sort \
        | awk '
            BEGIN { first = 1 }
            { printf "%s%s", first ? "" : ":", $0; first = 0 }
            END { print "" }
        '
}

matrix-save-shell-state() {
    if [ "${MATRIX_ENV_ACTIVE:-0}" = "1" ]; then
        return 0
    fi
    export MATRIX_SAVED_CWD="$PWD"
    export MATRIX_SAVED_PATH="${PATH:-}"
    export MATRIX_SAVED_PYTHONPATH="${PYTHONPATH:-}"
    export MATRIX_ENV_ACTIVE=1
}

matrix-ensure-gui-deps() {
    if [ "$MATRIX_AUTO_INSTALL_GUI_DEPS" = "0" ]; then
        return 0
    fi
    python - <<'PY' >/dev/null 2>&1
import matplotlib
import scipy
import PySide6
PY
    if [ $? -eq 0 ]; then
        return 0
    fi
    echo "Dipendenze MATRIX GUI/DVR mancanti: installo $MATRIX_GUI_DEPS..."
    matrix-pip-install-list "$MATRIX_GUI_DEPS"
}

matrix-pip-install-list() {
    _MATRIX_PIP_DEPS="$1" python - <<'PY'
import os
import shlex
import subprocess
import sys

deps = shlex.split(os.environ.get("_MATRIX_PIP_DEPS", ""))
if not deps:
    raise SystemExit(0)
raise SystemExit(subprocess.call([sys.executable, "-m", "pip", "install", *deps]))
PY
}

matrix-create-venv() {
    if [ "$MATRIX_AUTO_CREATE_VENV" = "0" ]; then
        return 1
    fi
    if [ -z "$MATRIX_VENV" ]; then
        echo "MATRIX_VENV non impostato: impossibile creare il virtualenv."
        return 1
    fi
    if ! command -v "$MATRIX_PYTHON" >/dev/null 2>&1; then
        echo "Python per MATRIX non trovato: $MATRIX_PYTHON"
        return 1
    fi
    echo "Creo virtualenv MATRIX: $MATRIX_VENV"
    mkdir -p "$(dirname "$MATRIX_VENV")" || return
    "$MATRIX_PYTHON" -m venv "$MATRIX_VENV"
}

matrix-activate-venv() {
    if [ -n "$MATRIX_VENV" ] && [ -f "$MATRIX_VENV/bin/activate" ]; then
        # shellcheck disable=SC1091
        source "$MATRIX_VENV/bin/activate"
        return 0
    fi
    if [ -f "$MATRIX_HOME/.venv/bin/activate" ]; then
        # shellcheck disable=SC1091
        source "$MATRIX_HOME/.venv/bin/activate"
        return 0
    fi
    if matrix-create-venv && [ -f "$MATRIX_VENV/bin/activate" ]; then
        # shellcheck disable=SC1091
        source "$MATRIX_VENV/bin/activate"
        return 0
    fi
    return 1
}

matrix-ensure-runtime-deps() {
    if [ "$MATRIX_AUTO_INSTALL_RUNTIME_DEPS" = "0" ]; then
        return 0
    fi
    python - <<'PY' >/dev/null 2>&1
import matplotlib
import numpy
import pandas
import pytest
import rdkit
import scipy
import sympy
PY
    if [ $? -eq 0 ]; then
        return 0
    fi
    echo "Dipendenze MATRIX runtime mancanti: installo $MATRIX_RUNTIME_DEPS..."
    python -m pip install --upgrade pip setuptools wheel || return
    matrix-pip-install-list "$MATRIX_RUNTIME_DEPS"
}

matrix-conda-hook() {
    if ! command -v conda >/dev/null 2>&1; then
        return 1
    fi
    local conda_base
    conda_base="$(conda info --base 2>/dev/null)" || return 1
    if [ -f "$conda_base/etc/profile.d/conda.sh" ]; then
        # shellcheck disable=SC1090
        source "$conda_base/etc/profile.d/conda.sh"
        return 0
    fi
    conda activate --help >/dev/null 2>&1
}

matrix-set() {
    matrix-save-shell-state
    if matrix-activate-venv; then
        :
    elif matrix-conda-hook && conda env list | awk '{print $1}' | grep -qx "$MATRIX_CONDA_ENV"; then
        conda activate "$MATRIX_CONDA_ENV"
    elif matrix-conda-hook && conda env list | awk '{print $1}' | grep -qx "matrix"; then
        conda activate matrix
    elif matrix-conda-hook && conda env list | awk '{print $1}' | grep -qx "oracle_26"; then
        conda activate oracle_26
    elif matrix-conda-hook && conda env list | awk '{print $1}' | grep -qx "oracle"; then
        conda activate oracle
    else
        echo "Warning: MATRIX environment not found, continuing with current Python."
    fi

    local package_path
    package_path="$(matrix-package-path)"
    if [ -n "$package_path" ]; then
        case ":${PYTHONPATH:-}:" in
            *":$package_path:"*) ;;
            *) export PYTHONPATH="$package_path${PYTHONPATH:+:$PYTHONPATH}" ;;
        esac
    fi

    export PATH="$MATRIX_HOME/tools:$MATRIX_HOME/scripts:$MATRIX_HOME/bin:$PATH"
    matrix-ensure-runtime-deps || return
    matrix-ensure-gui-deps || return
    cd "$MATRIX_HOME" || return
}

matrix-unset() {
    if [ -n "${VIRTUAL_ENV:-}" ]; then
        case "$VIRTUAL_ENV" in
            "$MATRIX_VENV"|"$MATRIX_HOME/.venv")
                if command -v deactivate >/dev/null 2>&1; then
                    deactivate
                fi
                ;;
        esac
    fi

    if [ -n "${CONDA_DEFAULT_ENV:-}" ] \
        && { [ "$CONDA_DEFAULT_ENV" = "$MATRIX_CONDA_ENV" ] || [ "$CONDA_DEFAULT_ENV" = "matrix" ]; } \
        && command -v conda >/dev/null 2>&1
    then
        conda deactivate
    fi

    if [ "${MATRIX_ENV_ACTIVE:-0}" = "1" ]; then
        export PATH="${MATRIX_SAVED_PATH:-$PATH}"
        if [ -n "${MATRIX_SAVED_PYTHONPATH+x}" ]; then
            if [ -n "$MATRIX_SAVED_PYTHONPATH" ]; then
                export PYTHONPATH="$MATRIX_SAVED_PYTHONPATH"
            else
                unset PYTHONPATH
            fi
        fi
        if [ -n "${MATRIX_SAVED_CWD:-}" ] && [ -d "$MATRIX_SAVED_CWD" ]; then
            cd "$MATRIX_SAVED_CWD" || return
        fi
    fi

    unset MATRIX_SAVED_CWD
    unset MATRIX_SAVED_PATH
    unset MATRIX_SAVED_PYTHONPATH
    unset MATRIX_ENV_ACTIVE
}

matrix-cli() {
    matrix-set || return
    python -m matrix "$@"
}

matrix-run() {
    matrix-set || return
    if python - <<'PY' >/dev/null 2>&1
import importlib.util
raise SystemExit(0 if importlib.util.find_spec("matrix_oracle.app") else 1)
PY
    then
        python -m matrix_oracle.app "$@"
    else
        python -m matrix "$@"
    fi
}

matrix-run-bg() {
    matrix-set || return
    local log="${1:-$MATRIX_HOME/logs/matrix_run.log}"
    mkdir -p "$(dirname "$log")"
    if python - <<'PY' >/dev/null 2>&1
import importlib.util
raise SystemExit(0 if importlib.util.find_spec("matrix_oracle.app") else 1)
PY
    then
        nohup python -m matrix_oracle.app > "$log" 2>&1 &
    else
        nohup python -m matrix > "$log" 2>&1 &
    fi
    echo "MATRIX avviato in background. PID: $! Log: $log"
}

matrix-test() {
    matrix-set || return
    if [ -z "$1" ]; then
        echo "Usage: matrix-test <test_script.py>"
        return 1
    fi
    python "$1"
}

matrix-test-all() {
    matrix-set || return
    python -m pytest "$MATRIX_HOME/tests" "$MATRIX_HOME/packages"
}

matrix-run-check() {
    matrix-set || return
    python - <<'PY'
import importlib

mods = ["matrix_core", "numpy", "scipy", "matplotlib", "pandas", "sympy", "PySide6", "rdkit"]
ok = True
for name in mods:
    try:
        mod = importlib.import_module(name)
        version = getattr(mod, "__version__", "n/a")
        print(f"[OK] {name} {version}")
    except Exception as exc:
        ok = False
        print(f"[FAIL] {name}: {exc}")
if ok:
    print("MATRIX runtime check: OK")
else:
    raise SystemExit(1)
PY
    local runtime_status=$?
    if [ "$runtime_status" -ne 0 ]; then
        return "$runtime_status"
    fi
    if command -v molden >/dev/null 2>&1; then
        echo "[OK] molden $(molden -h 2>/dev/null | awk 'NR==1 {print $1}')"
    else
        echo "[WARN] molden: not found in PATH"
    fi
    if [ "$(uname -s 2>/dev/null)" = "Darwin" ]; then
        if [ -d /Applications/Utilities/XQuartz.app ]; then
            echo "[OK] XQuartz installed"
        else
            echo "[WARN] XQuartz missing: Molden GUI windows require XQuartz on macOS"
        fi
    fi
}

matrix-install-gui-deps() {
    MATRIX_AUTO_INSTALL_GUI_DEPS=0 matrix-set || return
    matrix-pip-install-list "$MATRIX_GUI_DEPS"
}

matrix-install-runtime-deps() {
    MATRIX_AUTO_INSTALL_RUNTIME_DEPS=0 matrix-set || return
    python -m pip install --upgrade pip setuptools wheel || return
    matrix-pip-install-list "$MATRIX_RUNTIME_DEPS"
}

matrix-clean() {
    if [ "$1" = "outputs" ] || [ "$1" = "all" ]; then
        echo "Cleaning MATRIX output files outside tests/fixtures (.log .out .err)..."
        find . -path "./tests/fixtures" -prune -o -type f \
            \( -name "*.log" -o -name "*.out" -o -name "*.err" \) -delete
    else
        echo "Cleaning Python cache only..."
    fi

    find . -type d -name "__pycache__" -prune -exec rm -rf {} +
    find . -type f \( -name "*.pyc" -o -name "*.pyo" \) -delete

    echo "Done."
}

matrix-gic-corpus-list() {
    matrix-cli gicforge corpus --format paths --suffix .inp "$@"
}

matrix-gic-corpus-summary() {
    matrix-cli gicforge corpus "$@"
}

matrix-gic-corpus-audit() {
    matrix-cli gicforge corpus-audit "$@"
}

oracle-package-path() { matrix-package-path "$@"; }
oracle-save-shell-state() { matrix-save-shell-state "$@"; }
oracle-ensure-gui-deps() { matrix-ensure-gui-deps "$@"; }
oracle-pip-install-list() { matrix-pip-install-list "$@"; }
oracle-create-venv() { matrix-create-venv "$@"; }
oracle-activate-venv() { matrix-activate-venv "$@"; }
oracle-ensure-runtime-deps() { matrix-ensure-runtime-deps "$@"; }
oracle-conda-hook() { matrix-conda-hook "$@"; }
oracle-set() { matrix-set "$@"; }
oracle-unset() { matrix-unset "$@"; }
oracle-cli() { matrix-cli "$@"; }
oracle-run() { matrix-run "$@"; }
oracle-run-bg() { matrix-run-bg "$@"; }
oracle-test() { matrix-test "$@"; }
oracle-test-all() { matrix-test-all "$@"; }
oracle-run-check() { matrix-run-check "$@"; }
oracle-install-gui-deps() { matrix-install-gui-deps "$@"; }
oracle-install-runtime-deps() { matrix-install-runtime-deps "$@"; }
oracle-clean() { matrix-clean "$@"; }
oracle-gic-corpus-list() { matrix-gic-corpus-list "$@"; }
oracle-gic-corpus-summary() { matrix-gic-corpus-summary "$@"; }
oracle-gic-corpus-audit() { matrix-gic-corpus-audit "$@"; }

unset _matrix_env_file
unset _matrix_env_dir
