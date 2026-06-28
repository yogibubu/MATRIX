NH3 HF/cc-pVTZ nuclear quadrupole coupling golden fixture.

The snippets are intentionally minimal. They preserve the output rows needed to
test MATRIX quadrupole property promotion without committing full QM outputs.

Source runs:

- Gaussian GDV32, `#p hf/cc-pVTZ scf=tight prop=efg output=pickett nosymm`
- Molpro 2025.3, `{hf; expec,fg}`
- ORCA 6.1.1, `! HF cc-pVTZ TightSCF NoFrozenCore` with `%eprnmr fgrad`

The values are compared in the spectroscopic/Pickett convention. The transverse
`a/b` axes of NH3 are nearly degenerate, so cross-code tests compare the sorted
transverse pair plus the unique axis.
