# ADR 0002: Python And Strict Fortran77 Backends

Date: 2026-06-26

## Decision

Some MATRIX modules may intentionally provide both Python and strict Fortran77
implementations of the same scientific kernel.

This is allowed when the implementations are exposed as backends of one public
service, not as independent tools with divergent formats.

## Required Contract

Every dual-backend module must define:

- one public service API;
- one enriched XYZ section contract for reusable state;
- one set of normalized auxiliary input/output files when files are needed;
- one manifest schema;
- explicit backend selection in parameters;
- backend metadata in every run manifest;
- shared fixtures;
- Python/Fortran77 identity tests when exact identity is expected;
- tolerance-based regression tests when small numerical differences are
  expected.

## Why This Is Not Bad Duplication

The project benefits from dual implementations:

- Fortran77 preserves strict legacy behavior and validated numerical kernels.
- Python gives inspectable orchestration, testing, parsing and rapid evolution.
- Independent implementations are useful scientific cross-checks.
- Performance-critical kernels can remain in Fortran77 without forcing the rest
  of the system into backend-specific formats.

The duplication to avoid is not "same scientific formula in two languages"; it
is duplicate data models, duplicate parsers, duplicate section writers,
duplicate GUI logic and incompatible workflow contracts.

## Boundary

Allowed:

```text
matrix-neo service
  -> Python backend
  -> Fortran77 backend
  -> same #GIC / #SYCART sections
  -> same manifest keys
  -> same regression fixtures
```

Not allowed:

```text
GUI calls Fortran directly
script parses its own XYZ dialect
Python backend writes one schema
Fortran backend writes unrelated files with no adapter
tests cover only one backend
```

## Initial Dual-Backend Candidates

- GICForge coordinate construction and B-matrix evaluation.
- DVR path/grid solver kernels.
- GF/VPT2/VCI numerical kernels.
- Semiexperimental fitting kernels when legacy compatibility is required.
