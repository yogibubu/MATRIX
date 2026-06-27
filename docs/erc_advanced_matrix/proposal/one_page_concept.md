# One-Page ERC Concept Note Draft

Working title:

> MATRIX: knowledge-based molecular intelligence for predictive spectroscopy

## Core Claim

MATRIX will make predictive spectroscopy possible for light-element molecular
systems up to about 100 atoms, aiming at an accuracy that is currently
reachable only for systems of roughly 10 atoms. It will do so by representing
structure, topology, symmetry, internal coordinates, quantum data, force
fields, spectra and provenance as one coherent molecular state. The project
will replace fragile chains of specialised, disconnected calculations with a
physically interpretable framework where each observable is a controlled
projection of the same molecular object, implemented in a robust ecosystem that
can run from laptops to workstations and high-performance computing centres.

## Scientific Need

The frontier of molecular spectroscopy and nano-engineering is moving toward
larger, flexible, weakly bound and multifunctional systems. These are often too
large for brute-force high-level treatments and too chemically specific for
generic empirical models. Current workflows repeatedly parse and reinterpret
the same molecule, losing normal-mode identity, symmetry, topology, fragment
meaning and error provenance. This prevents reliable prediction and comparison
with experiment.

The decisive frontier is not to run the same protocol on ever larger systems.
For molecules around 1000 atoms, and for condensed phases, spectral resolution
and molecular specificity generally decrease, so the right tools are
multiscale: QM/MM, coarse-graining and continuum descriptions. MATRIX will
provide the accurate molecular physics and microsolvation-derived parameters
that those larger-scale models need.

The project also addresses a practical barrier to frontier spectroscopy:
advanced methods remain fragile when they require expert-only scripts,
machine-specific setups and manual transfers between codes, databases and
visualizers. MATRIX will expose the same scientific state through command-line,
graphical and web-aware interfaces, adapting execution to the available
hardware without changing the underlying model.

## MATRIX Idea

MATRIX will build a molecular reasoning framework based on an enriched,
versioned molecular state. This state will contain geometry, topology,
symmetry, synthons, fragments, generalized internal coordinates, B matrices,
force constants, Hessians, electronic data, rotational constants, spectra and
experimental links. NEO will provide the generalized internal-coordinate
language; the other MATRIX tools will use the same state rather than
duplicating parsers, topology analysis or physical assumptions.

The framework will communicate with external resources through documented
interfaces: experimental databases, molecular repositories, web services,
quantum-chemistry engines, visualization programs and remote computing
systems. The same calculation should be inspectable on a laptop, optimized on a
workstation and scaled on a computing centre.

MATRIX is not only an implementation and integration project. Several of its
tools are scientific contributions in their own right: localized and
symmetry-adapted generalized internal coordinates, chemically controlled
primitive-coordinate selection, electronic-structure PCS/composite strategies,
VPT2/VCI workflows tied to the same molecular state, and
local-mode/fragment-aware models for large light-element molecules. The
ecosystem exists to make these new methods usable and reproducible, not to
replace the theoretical advances.

## Accuracy Philosophy

MATRIX will not pursue accuracy by applying the same expensive level to every
term, nor by neglecting small contributions. Its strategy is to include every
physically meaningful contribution, compute each one at the level justified by
its weight in the final observable, and keep the result coherent and
upgradeable. This generalizes the logic of electronic-structure PCS/composite
schemes and three-level vibrational models: harmonic and anharmonic terms,
equilibrium structures and vibrational corrections, local and non-local
effects can be assigned different accuracies without losing consistency.

Intra-fragment interactions, which are often strongly correlated, and
inter-fragment interactions, which are usually weaker and more transferable,
will be treated at different levels not only for energies but also for
structures and force fields. Machine learning enters only at the end, as a
DeltaML correction for the residual random-like error that remains after this
knowledge-based decomposition.

## DeltaML Philosophy

Machine learning will be used only as DeltaML. MATRIX will first remove
systematic errors by knowledge-based models: improved internal coordinates,
symmetry adaptation, topology, fragment descriptors, physically meaningful
force constants, electrostatic and dispersion corrections, anharmonic models,
multi-level quantum data and intra-/inter-fragment accuracy partitioning. Only
the remaining residual, effectively random error will be learned. The goal is
not black-box prediction, but interpretable error completion after chemical
understanding has done all it can.

## Target Domain

The project will focus on molecules up to about 100 atoms, especially organic,
bio-inspired, atmospheric, supramolecular and nano-engineered systems. It will
avoid metals and heavy elements, allowing at most one heavy atom when
scientifically essential. This boundary targets sustainable, abundant-element
molecular science and keeps the dominant physics in a regime where
spectroscopic interpretation remains possible.

## Breakthrough Objectives

1. Define a persistent molecular state that preserves identity across
   structure, spectroscopy, thermochemistry and electronic observables.
2. Develop symmetry-, topology- and fragment-aware generalized internal
   coordinates with analytic B matrices and robust reduction protocols.
3. Create new localized/symmetry-adapted GIC, electronic-structure
   PCS/composite, VPT2/VCI and local-mode methodologies as scientific
   advances, not only as software modules.
4. Generalize PCS/three-level accuracy decomposition so that all contributions
   are included coherently and can be upgraded independently.
5. Build predictive vibrational, rotational, electronic and thermochemical
   workflows for systems approaching 100 atoms.
6. Create knowledge-based error decomposition and DeltaML residual correction.
7. Transfer accurate molecular information to QM/MM, coarse-grained and
   continuum models for larger systems and condensed phases.
8. Build a robust, user-friendly and hardware-aware ecosystem with command-line,
   GUI, web/database and high-performance-computing interfaces.
9. Validate MATRIX against demanding corpora and experimental databases,
   producing publication-ready spectra and traceable predictions.

## High-Risk / High-Gain Point

The high-risk question is whether molecular information can be preserved and
propagated rigorously enough that large-system spectroscopy becomes predictive
without losing interpretability. The gain would be a new computational paradigm:
not a collection of codes, but a chemically aware molecular intelligence
framework for sustainable molecular design and high-resolution spectroscopy.
