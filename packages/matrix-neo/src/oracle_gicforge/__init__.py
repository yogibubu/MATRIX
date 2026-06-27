"""Compatibility alias for :mod:`matrix_neo`."""

from importlib import import_module
import sys

_module = import_module("matrix_neo")
globals().update(_module.__dict__)
sys.modules[__name__] = _module
