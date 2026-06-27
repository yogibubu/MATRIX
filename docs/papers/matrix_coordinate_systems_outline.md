# MATRIX Automatic Coordinate Systems: Paper Outline

## Working Title

Automatic, symmetry-adapted and locality-preserving internal coordinates for
molecular spectroscopy and quantum-chemical workflows.

## Core Thesis

MATRIX provides an automatic coordinate construction strategy that combines four
requirements that are usually handled separately:

- locality of the chemical description;
- coordinate-type separation during reduction and symmetrization;
- exact point-group symmetry adaptation when symmetry is present;
- protected special coordinates for rings, fragments, centers and weak contacts.

The goal is not only to generate non-redundant internal coordinates, but to
generate chemically meaningful coordinates that remain suitable for geometry
optimization, least-squares fitting, GF analysis and force-field development.

## Scientific Novelty

- Localized GICs are generated from a single saved topology and reused by all
  MATRIX tools.
- Ring puckering, fused-ring torsions and butterfly coordinates are treated as
  protected coordinate families before non-redundant reduction.
- Symmetry adaptation is performed by point-group operations without mixing
  coordinate types.
- Non-covalent complexes can be described by two complementary models:
  fragment-center coordinates with relative translations/rotations, or
  pseudo-bond/H-bond coordinates added to the internal-coordinate graph.
- Special coordinates remain compatible with analytic Wilson B-matrix
  construction.

## Coordinate Models

### Fragment Coordinate Model

Each fragment keeps its own internal coordinates. Inter-fragment motion is
described by relative center distances, translations and orientations. This is
appropriate for weakly bound complexes, host-guest systems and fragment-based
fits where rigid-body motions are chemically meaningful.

### Pseudo-Bond / H-Bond Model

Weak contacts are represented by pseudo-bonds. The graph is connected through
the closest inter-fragment atom pair, or through detected X-H...Y hydrogen bonds
when chemically available. Once the pseudo-bond is added, standard internal
coordinates are generated and reduced as for an ordinary molecule. Pseudo-bonds
do not create artificial rings.

## Required Numerical Demonstrations

- Small non-covalent probes: formic-acid water, water dimer and mixed H-bonded
  complexes.
- Ring systems: benzene, pyrrole, azulene, pyrene and fused/polycyclic cases.
- Bridged/spiro systems: norbornane, norbornene, norbornadiene, thujone,
  ribose, cubane and spiro probes.
- Metal/center cases: ferrocene eclipsed and staggered, including D5h/D5d
  symmetry closure.
- Python/Fortran consistency checks on primitive sets, GICs and B rows.

## Key Comparisons

- Standard redundant internal coordinates versus MATRIX reduced GICs.
- Fragment-coordinate model versus pseudo-bond model for the same weak complex.
- Unsymmetrized versus point-group symmetrized GICs.
- Python and Fortran implementations on the same regression corpus.

## Expected Outputs

- Stable GIC definitions and Gaussian ReadGIC blocks.
- Analytic B matrices for all supported coordinate families.
- Total-symmetric coordinate subsets for optimization and least-squares fitting.
- Reproducible xyzin sections documenting topology, symmetry, fragment mode,
  pseudo-bonds and frozen GICs.

## Open Technical Points

- Decide when the automatic default should prefer fragment coordinates or
  pseudo-bonds for weak complexes.
- Extend H-bond selection beyond closest-pair fallback using the Merlino/BDPCS3
  donor-H-acceptor logic.
- Define publication-quality diagnostics for coordinate locality, conditioning
  and symmetry purity.
- Complete periodic Python-versus-Fortran B-row audits for special coordinates.
