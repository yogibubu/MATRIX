# Scale And Accuracy Strategy

## Scientific Positioning

MATRIX is not a "make everything bigger" project. Its central claim is to move
the boundary of high-resolution, interpretable molecular spectroscopy from
systems of roughly 10 atoms to light-element systems approaching 100 atoms.

This size range is large enough for realistic molecular function,
supramolecular interactions, nano-engineering motifs, conformational
flexibility and weak interactions. It is still small enough that molecular
identity, normal-mode assignment, symmetry, topology and internal-coordinate
interpretation remain scientifically meaningful.

## Why Not 1000 Atoms

For systems around 1000 atoms, the experimental and theoretical problem changes
character. Spectral resolution usually decreases, individual molecular
features become less separable, environmental broadening becomes dominant and
the full high-resolution isolated-molecule description is often no longer the
right level of theory.

MATRIX should therefore not claim that the same spectroscopic protocol applies
unchanged to all scales. Instead, the project should define a scale hierarchy:

- up to about 10 atoms: current high-accuracy reference regime;
- up to about 100 atoms: MATRIX frontier regime for high-resolution,
  knowledge-based predictive spectroscopy;
- around 1000 atoms and condensed phases: multiscale regime using QM/MM,
  coarse-grained and continuum descriptions.

## Role In Multiscale Models

MATRIX remains relevant beyond 100 atoms by providing accurate molecular
building blocks:

- force constants and internal-coordinate parameters;
- physically corrected charges, polarizabilities and dispersion parameters;
- local-mode and fragment descriptors;
- conformer and microsolvation data;
- benchmark spectra and uncertainty estimates.

For condensed phases, the proposal should emphasize accurate microsolvation as
the bridge. Local solvent clusters can be treated with MATRIX-level accuracy to
derive MM parameters and environmental corrections that feed larger QM/MM,
coarse-grained or continuum simulations.

## ERC Message

The scale boundary is a strength, not a limitation. MATRIX targets the regime
where high-resolution spectroscopy is currently blocked by loss of molecular
information, not the regime where the physics itself demands lower-resolution
statistical or continuum descriptions.
