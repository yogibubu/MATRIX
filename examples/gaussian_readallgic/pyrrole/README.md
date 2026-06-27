# Pyrrole Gaussian ReadAllGIC Example

This directory contains a successful Gaussian DV J32 run generated from MATRIX
NEO/GICForge symmetrized GICs.

Files:

- `pyrrole.gjf`: Gaussian input using `opt=readallgic freq`.
- `pyrrole.log`: Gaussian output showing successful optimization and frequency
  calculation.

The example demonstrates that:

- Gaussian reads the MATRIX-generated ReadAllGIC block.
- The full 24-coordinate nonredundant basis is present (`NTRed=24`,
  `NRank=24`).
- Only the totally symmetric `A1` coordinates remain active.
- Non-totally-symmetric coordinates, including `B2` stretches/bends,
  `A2/B1` ring puckering and `A2/B1` out-of-plane coordinates, are marked
  `Frozen`.
- Gaussian preserves `C2V` symmetry, finds a stationary point and reports
  `NImag=0`.

The example is covered by
`tests/test_matrix_gaussian_readallgic_examples.py`, so CI can check the
Gaussian-side ReadAllGIC contract without launching Gaussian.

Rebuild the MATRIX sections from the optimized Gaussian log and run GF/PED with:

```bash
matrix link preprocess pyrrole.log pyrrole.xyzin --source-kind gaussian
matrix validate pyrrole.xyzin
matrix neo build pyrrole.xyzin --symmetrize
matrix gaussian promote-log-hessian pyrrole.log pyrrole.xyzin
matrix gf --xyzin pyrrole.xyzin --symmetry-blocks
```

The `.gjf` ReadAllGIC block is not imported as the MATRIX GIC definition.
NEO/MATRIX rebuilds `#GIC` from the final optimized log geometry, molecular
topology and symmetry policy. The Gaussian log is then also used as an external
QM result: `promote-log-hessian` writes `#CARTESIAN_HESSIAN` and, when the
printed table is present, `#NORMAL_MODES`.

For Gaussian frequency jobs the printed normal coordinates are emitted in the
Gaussian normal-coordinate frame, while the Cartesian Hessian is restored to
the archive/original axes. The adapter stores the archive/original geometry and
rotates the printed normal modes to that frame when the geometry fit is
unambiguous. GF/PED reports both a geometry check and a frequency check against
the source Gaussian frequencies; the GF frequencies reproduce the Gaussian
frequencies within the rounding of the printed log matrix.

If `--symmetry-blocks` is requested, GF/PED first verifies that the internal
`F` and `G` matrices are block diagonal in the stored irrep labels. If
off-block couplings are larger than numerical roundoff, GF/PED stops because
correctly symmetrized coordinates cannot couple different irreps. The pyrrole
fixture exercises the strict block solve and checks that every PED column is
finite, non-negative, normalized to 100%, and confined to its symmetry block.

Use `matrix gaussian promote-log-hessian --no-normal-modes pyrrole.log
pyrrole.xyzin` when only the Cartesian Hessian section should be refreshed.
