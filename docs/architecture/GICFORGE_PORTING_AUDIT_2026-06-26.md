# GICForge Porting Audit

Date: 2026-06-26

## Rule

ORACLE ports Merlino3.0 first. A behavior can change only when Merlino3.0 has a
documented fragility or bug, and the ORACLE change must be covered by a
regression test and noted in the method documentation.

## Sources Reviewed

- `/Users/vincenzobarone/merlino3.0/merlino_gic/model.py`
- `/Users/vincenzobarone/merlino3.0/merlino_gic/gic_symmetry.py`
- `/Users/vincenzobarone/merlino3.0/merlino_fit/survibfit/primitives.py`
- `/Users/vincenzobarone/merlino3.0/merlino_fit/survibfit/transforms.py`
- `/Users/vincenzobarone/merlino3.0/fortran/gicforge`
- `/Users/vincenzobarone/Desktop/newmsr_overleaf/gicforge_manual.tex`
- `/Users/vincenzobarone/Desktop/newmsr_overleaf/main.tex`
- `/Users/vincenzobarone/Desktop/newmsr_overleaf/manuals_README.md`

## Porting Matrix

| Area | Merlino3.0 source | ORACLE status | Action |
| --- | --- | --- | --- |
| Ordinary primitives | `survibfit/primitives.py`, `mkprim.f` | Present: stretch, bend, linear bend, torsion, out-of-plane | Keep parity tests against analytic/FD B rows. |
| Ring bend/torsion blocks | `survibfit/transforms.py`, `mkcyc.f`, `mksalc.f`, Overleaf GICForge manual | Port started: `CYCLIC_BEND`, `CYCLIC_TORSION`, `CONDENSED_RING_TORSION`, `BUTTERFLY` families now exist in Python and remain separate for symmetry/reduction | Next: implement full Merlino ring SALC/puckering combinations `QPck/PhiP/RPck` instead of only family classification. |
| Butterfly coordinates | `mksalc.f:BtFly`, `survibfit/transforms.py:ring_butterfly_u` | Port started: fused-ring shared-bond torsions are classified as `BUTTERFLY` | Next: compare coefficient vectors with Merlino for polycyclic corpus fixtures. |
| Non-redundant reduction | `gicprune.f`, `locsvd.f`, `transforms.py` | Present: analytic B-row MGS, protected special-first policy, ring families not mixed with ordinary blocks | Next: add corpus-level count/golden tests for ring systems. |
| Symmetrization | `gic_symmetry.py`, `symmetry_global.py`, `gic_type_symmetry.f`, `symm.f` | Present: Merlino label-only parity plus matrix projectors; no type mixing; total-symmetric subset stored | Next: full strict Fortran projector diagnostics, not only local SALC/static source checks. |
| Fragment/TRIC coordinates | Merlino TRIC roadmap, Overleaf manual, geomeTRIC reference model | Present: fragment center distance, fragment center-atom distance, translations, orientations, analytic B rows, Gaussian symbolic export | Next: add atom-frame angle/torsion and center-frame tilt/orientation modes. |
| Ring/bond/interaction centers | Fragment roadmap, topology ring docs | Present: bond centers, ring centers, atom-center distance candidates, analytic chain-rule B row | Next: add center-angle, center-torsion and hapticity/coordination center scoring. |
| Python/Fortran parity | `doc/GIC_PYTHON_FORTRAN_COMPARISON.md`, Fortran GICForge sources | Partial: legacy source vendored and compiled; `frag_tric_bmat.f` mirrors special-coordinate B rows; tests assert legacy ring/symmetry sources remain present | Next: executable parity harness comparing primitive lists, GIC labels and B rows on corpus cases. |
| `xyzin` standalone restart | Merlino frozen schema docs | Present: frozen `#GIC`, `#SYCART`, report, B-matrix from saved sections | Next: extend downstream GF/SEFit/Thermo/VPT2 commands to consume the saved GIC section directly. |

## New Regression Anchors

- `tests/test_oracle_gicforge.py::test_gicforge_label_only_characters_match_merlino3`
  freezes Merlino3.0 label-only character behavior.
- `tests/test_oracle_gicforge.py::test_gicforge_classifies_ring_and_butterfly_primitives_like_merlino`
  freezes the Merlino ring-family separation for cyclic bend/torsion,
  condensed-ring torsion and butterfly torsion.
- `tests/test_oracle_fortran_gicforge.py::test_legacy_merlino_ring_and_butterfly_blocks_remain_reference`
  prevents accidental removal of the strict Fortran77 ring/butterfly reference
  routines.

## Known Merlino Fragilities To Correct In ORACLE

- Label-only symmetry depends on operation spelling. ORACLE preserves it for
  compatibility, but uses matrix-classified projectors when operation matrices
  are available.
- Ring primitive generation in Merlino is split across Python and Fortran paths.
  ORACLE must centralize the saved topology/ring source in `xyzin` and make all
  backends consume the same ring records.
- Fragment/center derivatives were not one shared service in Merlino. ORACLE
  centralizes these rows in `oracle-gicforge` and mirrors them in
  `frag_tric_bmat.f`.

## Immediate Remaining Work

1. Port the full ring-puckering coordinate construction (`QPck`, `PhiP`, `RPck`)
   from the Merlino/Overleaf specification.
2. Add golden corpus cases for pyridine, pyrene, naphthalene/coronene-like fused
   rings and norcamphor/testosterone-like bridged systems.
3. Extend the strict Fortran executable harness to emit comparable primitive
   lists and B rows for the same frozen `#GIC` contract.
4. Add center-angle, center-torsion and atom-frame coordinates for
   ring/metal/H-bond cases after the topology-center records are finalized.
