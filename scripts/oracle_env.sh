# Compatibility wrapper for legacy ORACLE shell helpers.
#
# New shells should source:
#   source /Users/vincenzobarone/MATRIX/scripts/matrix_env.sh
#
# This file remains so existing .bashrc/.zshrc entries keep working during the
# MATRIX transition.

if [ -n "${BASH_SOURCE:-}" ]; then
    _oracle_env_wrapper_file="${BASH_SOURCE[0]}"
elif [ -n "${ZSH_VERSION:-}" ]; then
    _oracle_env_wrapper_file="${(%):-%x}"
else
    _oracle_env_wrapper_file="$0"
fi

_oracle_env_wrapper_dir="$(cd "$(dirname "$_oracle_env_wrapper_file")" 2>/dev/null && pwd)"

# shellcheck disable=SC1091
source "$_oracle_env_wrapper_dir/matrix_env.sh"

unset _oracle_env_wrapper_file
unset _oracle_env_wrapper_dir
