# MATRIX Legacy Repository Audit, 2026-06-28

Scope: `/Users/vincenzobarone/merlino3.0`, `/Users/vincenzobarone/merlino_fit`
and `/Users/vincenzobarone/newmerlin`.  `merlino3.0` remains the reference
source and should not be deleted.

## `/Users/vincenzobarone/merlino_fit`

This is not primarily a MATRIX software repository.  It contains CV/manuscript
material, PNO/LPNO paper material, CeDiTT4 material, an Angewandte H-bond
manuscript work area and a small diagnostic script.

Useful for MATRIX:

- `scripts/gaussian_one_mode_probe.py`;
- `scripts/README_gaussian_one_mode_probe.md`;
- `alpha_resonances_data/formic_rDSD.log`;
- `alpha_resonances_data/formic.fchk`;
- `alpha_resonances_data/malon_parent.log`;
- `alpha_resonances_data/malon_parent.fchk`.

Reason: the script is a validated Gaussian-consistent one-mode probe for
normal-coordinate alignment, `dI/dQ`, `I(Q)`, `1/I(Q)` and a 1D
variational/perturbative comparison.  MATRIX has VPT2/VCI and DVR machinery,
but this exact rovibrational one-mode diagnostic is not yet a public MATRIX
workflow.

Recommendation: do not treat `merlino_fit` as a MATRIX dependency.  Archive or
move it as a manuscript/work folder.  Before deleting, either keep the whole
folder elsewhere or copy the one-mode probe and its two small Gaussian
log/fchk pairs into a future MATRIX examples or research-notes area.

## `/Users/vincenzobarone/newmerlin`

This folder is a mixed historical working tree: old Fortran sources, Python
SMILES experiments, many Gaussian/Merlino inputs and very large calculation
artifacts.  It is not a git repository.  Its size is dominated by binary
checkpoint files, especially `phenylpeniciline/phenylpeniciline.chk`
(`~22 GB`).

Useful material not fully covered by current MATRIX:

- historical Fortran one-dimensional variational code: `vci1d.f` and
  `allvib.f`;
- kinetics/rate prototypes: `rates.f`;
- substitution/Kraitchman prototype: `source/newrot/newkra.f`;
- internal-coordinate derivative experiments: `der12.F`;
- historical SMILES/Z-matrix experiments: `smiles_to_zmat.f`, `newsmi.f`,
  `setZ.f`;
- old notes in `source/merlino_notes` and `source/doc/merlino_notes`;
- regression inputs missing from the current MATRIX corpus:
  `acetylene_smile.inp`, `azulene_smile.inp`, `cyclopentadiene_msr.inp`,
  `glycine.inp`, `glycine_BDPCS3.inp`, `glycine_HPCS2.inp`,
  `glycine_new.inp`, `norbornadiene_PCS0.inp`, `norbornadiene_PCS1.inp`,
  `testmsr.inp`, `tiopronin_smile_g16.inp`, `tiopronin_smile_gdv.inp`.

Material already covered better elsewhere:

- most GIC/ring/symmetry ideas are superseded by `merlino3.0` plus MATRIX
  NEO/GICForge and the current golden corpus;
- SMILES import should continue through MATRIX LINK/RDKit, not the old
  Fortran SMILES prototype;
- many `.chk`, `.log`, `.out`, `.fchk` files are calculation artifacts rather
  than source.

Recommendation: do not delete `newmerlin` blindly.  First extract the source
files, notes and missing `.inp` fixtures listed above into a small archive or
MATRIX TODO/import branch.  After that extraction, the bulky calculation
artifacts and duplicated working directories can be deleted or moved to cold
storage.  If disk space is urgent, the large `.chk` files are the first safe
cleanup candidates, provided the corresponding calculations do not need to be
restarted.

## Bottom Line

- Keep `merlino3.0`.
- `merlino_fit`: not needed as a MATRIX repo after preserving manuscript
  material and the one-mode probe.
- `newmerlin`: keep only until the listed source snippets and missing fixtures
  have been extracted; then it can be archived or cleaned aggressively.
