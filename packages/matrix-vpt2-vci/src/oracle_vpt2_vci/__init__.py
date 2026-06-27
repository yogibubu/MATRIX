"""Compatibility alias for :mod:`matrix_vpt2_vci`."""

from importlib import import_module
import sys

_module = import_module("matrix_vpt2_vci")
globals().update(_module.__dict__)
sys.modules[__name__] = _module
