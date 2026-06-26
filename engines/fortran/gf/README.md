# ORACLE GF Fortran Backend

This directory vendors the Merlino3.0 fixed-form Fortran77 Wilson-GF helper.

- `legacy_merlino/gf_core.f` is copied from
  `/Users/vincenzobarone/merlino3.0/fortran/vpt2_vci/gf_core.f`.
- The Python `oracle_gf.harmonic.solve_wilson_gf` path mirrors the same
  symmetric Wilson-GF construction.

The source is intentionally kept close to Merlino3.0. Change it only when a
specific Merlino fragility is identified and the Python path is updated in the
same way.
