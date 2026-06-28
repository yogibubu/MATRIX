# MATRIX legacy source staging

This directory keeps selected legacy material that should not be lost while
MATRIX replaces the older development trees.  Files here are not active runtime
code unless they are explicitly ported into a package under `packages/`.

## `merlino_fit_one_mode_probe`

Source: `/Users/vincenzobarone/merlino_fit/scripts/`

- `gaussian_one_mode_probe.py`
- `README_gaussian_one_mode_probe.md`

These are useful for future anharmonic/local-mode probes.  The large data files
under `/Users/vincenzobarone/merlino_fit/alpha_resonances_data/` were not copied
into the repository:

- `formic_rDSD.log`
- `formic.fchk`
- `malon_parent.log`
- `malon_parent.fchk`

They remain good candidate regression data, but they should be imported only if
we decide that the extra repository size is justified.

## `newmerlin`

Source: `/Users/vincenzobarone/newmerlin/`

The staged Fortran files cover material not yet fully absorbed by MATRIX:

- `vci1d.f`, `allvib.f`, `rates.f`: legacy anharmonic/rate prototypes.
- `newkra.f`: Kraitchman-related rotational material.
- `der12.F`: derivative/probing utility.
- `smiles_to_zmat.f`, `newsmi.f`, `setZ.f`: old SMILES/Z-matrix experiments.

The staged notes capture implementation reminders from `source/merlino_notes`
and `source/doc/merlino_notes`.

Small input files from `newmerlin/tests` and `newmerlin/source/tests` are staged
under `tests/fixtures/test_molecules/newmerlin_legacy/`.  They can be promoted
into official golden corpora once the relevant workflows are validated.
