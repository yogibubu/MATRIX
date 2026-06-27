"""Source-tree launcher shim for :mod:`oracle_gui`."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "packages" / "oracle-gui" / "src"
PKG_ROOT = SRC_ROOT / "oracle_gui"

for src in sorted((ROOT / "packages").glob("*/src")):
    text = str(src)
    if text not in sys.path:
        sys.path.insert(0, text)
if str(PKG_ROOT) not in __path__:
    __path__.insert(0, str(PKG_ROOT))

_real_init = PKG_ROOT / "__init__.py"
exec(compile(_real_init.read_text(encoding="utf-8"), str(_real_init), "exec"), globals())
