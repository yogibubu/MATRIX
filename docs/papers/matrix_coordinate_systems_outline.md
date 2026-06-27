# MATRIX Automatic Coordinate Systems: Paper Draft

## Working Title

Automatic, symmetry-adapted and locality-preserving internal coordinates for
molecular spectroscopy and quantum-chemical workflows.

## Draft Abstract

We present an automatic construction of generalized internal coordinates for
spectroscopic and quantum-chemical workflows.  The method starts from a single
validated molecular state containing geometry, topology, symmetry and optional
fragment information, and produces frozen coordinate definitions together with
analytic Wilson B matrices.  The construction is designed to preserve chemical
locality while avoiding the loss of physically important coordinates during
non-redundant reduction.  Coordinates are generated and reduced in homogeneous
families: stretches, bends, torsions, out-of-plane/improper coordinates, ring
puckering coordinates, butterfly coordinates and protected special coordinates.
Point-group symmetry adaptation is carried out by projector operations without
mixing different coordinate types.  Weak complexes can be described either by
protected fragment-center translation/rotation coordinates or by hydrogen-bond
pseudo-bonds that connect fragments and then use ordinary internal coordinates.
The validity of the symmetry assignment is tested downstream: in a
symmetry-adapted GIC basis, Wilson \(G\) and transformed force-constant \(F\)
matrices must be block diagonal by irreducible representation, and potential
energy distributions must remain within the same symmetry block.
The resulting framework unifies coordinate construction, symmetry adaptation
and reusable B-matrix evaluation for geometry optimization, least-squares
structural refinement, GF analysis and force-field development.

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
- Ring-puckering coordinates are validated in the physical B-row source space
  when the selected `RPck` source subspace is lower-dimensional than the
  redundant endocyclic torsion set. This keeps the symmetry projector tied to
  the same derivatives used by GF, fitting and optimization.
- Symmetry labels are checked by an independent GF invariant: cross-irrep
  couplings in \(F\) or \(G\) are contract failures, not a downstream
  diagonalization choice.
- Weak non-covalent contacts are handled without contaminating the covalent
  topology: pseudo-bonds are local graph edges for coordinate generation, not
  new covalent bonds.
- The same molecule can be represented by alternative coordinate models whose
  B ranks are tested against the same vibrational space.

## Algorithmic Contract

1. Read or construct a validated `xyzin` molecular state.
2. Build the covalent topology once and save it.
3. Determine point-group symmetry operations and atom permutations.
4. Generate primitive coordinates by coordinate family.
5. Add protected special-coordinate families when requested by fragments,
   virtual centers, rings or metal/center contacts.
6. Reduce only within physically meaningful classes using analytic B-matrix
   rank tests.
7. Apply point-group projector symmetrization without mixing coordinate types.
8. Save frozen GICs, primitive coefficients, symmetry labels, total-symmetric
   subsets and B-matrix metadata.
9. Validate the frozen coordinate model by rank, family-count, symmetry-label,
   B-row and downstream GF/PED invariants.

This contract is intended to make downstream tools deterministic: GF, SEFit,
anharmonic modules and GUI workflows consume the frozen coordinate schema
instead of rebuilding their own coordinate parser or primitive set.

## Coordinate Models

### Fragment Coordinate Model

Each fragment keeps its own internal coordinates. Inter-fragment motion is
described by relative center distances, translations and orientations. This is
appropriate for weakly bound complexes, host-guest systems and fragment-based
fits where rigid-body motions are chemically meaningful.

### Pseudo-Bond / H-Bond Model

Weak contacts are represented by pseudo-bonds. The graph is connected through
detected X-H...Y hydrogen bonds when chemically available, using the
Merlino/BDPCS3 donor, acceptor, distance and angular criteria. If no chemically
valid hydrogen bond connects fragments, the closest inter-fragment atom pair is
used as a fallback. Once the pseudo-bond is added, standard internal coordinates
are generated and reduced as for an ordinary molecule. Pseudo-bonds do not
create artificial rings.

For the formic-acid--water probe, the fragment model gives protected fragment
distance/translation/orientation coordinates, while the pseudo-bond model
detects the O-H...O contact and builds only ordinary stretches, bends and
torsions. Both models reach the full rank of the same 18-dimensional
vibrational space.

## Proposed Figures

- Workflow diagram: `xyzin` state -> topology/symmetry -> primitive families ->
  reduction -> projector symmetrization -> frozen GIC/B matrix.
- Coordinate-family diagram showing separate reduction blocks and forbidden
  type mixing.
- Formic-acid--water comparison: fragment-coordinate model versus H-bond
  pseudo-bond model.
- Ring-coordinate diagram for isolated and fused rings, including RPck/QPck/PhiP
  and butterfly coordinates.
- Ferrocene example showing center/ring coordinates and D5h/D5d closure.

## Proposed Tables

- Coordinate families, reduction class, Gaussian expression and analytic
  B-matrix status.
- Non-covalent model comparison for weak complexes.
- Regression corpus: molecule, point group, rank, special coordinate families,
  Python/Fortran status.
- Symmetry projector validation: group, operation count, closed permutations,
  total-symmetric coordinates.
- Cross-module invariant checks: final log geometry, Hessian source, GF
  symmetry blocks, off-block \(F/G\) maximum, PED normalization and irrep
  purity.

## Required Numerical Demonstrations

- Small non-covalent probes: formic-acid water, water dimer and mixed H-bonded
  complexes.
- Ring systems: benzene, pyrrole, azulene, pyrene and fused/polycyclic cases.
- Bridged/spiro systems: norbornane, norbornene, norbornadiene, thujone,
  ribose, cubane and spiro probes.
- Metal/center cases: ferrocene eclipsed and staggered, including D5h/D5d
  symmetry closure.
- Python/Fortran consistency checks on primitive sets, GICs and B rows.
- Gaussian ReadAllGIC closed-loop checks: MATRIX-generated symmetrized GICs,
  Gaussian optimization/frequency output, final log geometry recovery,
  Hessian promotion, strict symmetry-block GF/PED and symmetry-pure PED.

## Key Comparisons

- Standard redundant internal coordinates versus MATRIX reduced GICs.
- Fragment-coordinate model versus pseudo-bond model for the same weak complex.
- Unsymmetrized versus point-group symmetrized GICs.
- Python and Fortran implementations on the same regression corpus.
- Fragment model versus pseudo-bond model for the same non-covalent complex:
  same target rank and B rank, different physical coordinate interpretation.

## Expected Outputs

- Stable GIC definitions and Gaussian ReadGIC blocks.
- Analytic B matrices for all supported coordinate families.
- Total-symmetric coordinate subsets for optimization and least-squares fitting.
- Reproducible xyzin sections documenting topology, symmetry, fragment mode,
  pseudo-bonds and frozen GICs.
- Failure diagnostics that identify the violated invariant: rank, family count,
  operation closure, B-row independence, cross-irrep \(F/G\) coupling or PED
  symmetry leakage.

## Open Technical Points

- Decide when the automatic default should prefer fragment coordinates or
  pseudo-bonds for weak complexes.
- Broaden the H-bond detector beyond N/O acceptors when chemically justified,
  while keeping the Merlino/BDPCS3 default reproducible.
- Define publication-quality diagnostics for coordinate locality, conditioning
  and symmetry purity.
- Complete periodic Python-versus-Fortran B-row audits for special coordinates.

## Current MATRIX Regression Anchors

- Benzene: D6h symmetry, ring-puckering projector behavior.
- Pyrrole: C2v ring coordinates, non-total-symmetric RPck/out-of-plane labels,
  Gaussian ReadAllGIC closed-loop GF/PED with final log geometry and strict
  cross-irrep \(F/G\) rejection.
- Azulene: nonbenzenoid/fused-ring behavior.
- Norbornane/norbornene/norbornadiene/norcamphor/thujone/ribose/cubane:
  bridged, cyclic and high-connectivity cases.
- Ferrocene eclipsed and staggered: D5h/D5d closure and metal-to-ring-center
  coordinates.
- Formic-acid--water: fragment coordinates and H-bond pseudo-bond coordinates
  on the same non-covalent complex.

## Near-Term Writing Plan

1. Freeze terminology: GIC, primitive, protected special coordinate, fragment
   mode, pseudo-bond, projector block.
2. Extract one clean numerical table from the regression suite.
3. Generate one figure for the formic-acid--water dual coordinate model.
4. Generate one figure for ring/fused-ring coordinate families.
5. Add a Python/Fortran parity subsection explaining what is identical now and
   what remains backend-specific.
