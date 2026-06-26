# ORACLE GICForge Fortran Backend

This directory contains the strict Fortran77-side building blocks used while
porting the legacy Merlino GICForge backend into ORACLE.

`frag_tric_bmat.f` is the first shared fragment-coordinate library. It mirrors
the native Python GICForge definitions for:

- `FC_DIST`: fragment-center / fragment-center distance;
- `FCA_DIST`: fragment-center / atom distance;
- `FTRANS`: Cartesian component of a fragment-center displacement;
- `FROT`: TRIC/geomeTRIC-style exponential-map rotation component between two
  fragment local frames.

The `FROT` B row is analytic: derivatives are propagated through the centroid,
local frame, relative rotation matrix, quaternion and exponential-map
small-rotation limit. It is intended to be called from the imported legacy
`MkBNew` path instead of adding another local derivative implementation.
