# Part B1 - Extended Synopsis Draft

Working title:

> MATRIX: a unified molecular-intelligence framework for predictive
> spectroscopy, structure and dynamics

Target length: follow the official ERC-2026-AdG template.

## One-Sentence ERC Claim

MATRIX will transform molecular spectroscopy from a collection of disconnected
specialised calculations into a unified, symmetry-aware and uncertainty-aware
scientific framework that can predict, interpret and design molecular
observables for light-element systems up to about 100 atoms with an accuracy
now reachable only for much smaller molecules, across structure, vibrations,
rotations, electronic states and reactivity, in an ecosystem that adapts from
laptops to workstations and high-performance computing centres.

## Frontier Problem

Current molecular modelling workflows often lose information when moving from
geometry to topology, internal coordinates, force fields, spectra,
thermochemistry and electronic observables. The same molecule is repeatedly
parsed, reinterpreted and approximated by independent tools. This fragmentation
limits reliability for large, fluxional, weakly bound and highly anharmonic
systems, exactly where modern spectroscopy and molecular design need predictive
power.

The practical frontier is not the smallest benchmark molecules, but chemically
realistic systems approaching 100 atoms: conformationally rich organic,
bio-inspired, atmospheric, supramolecular and nano-engineered molecules for
which high-resolution spectroscopy is still possible but current fragmented
workflows become fragile. MATRIX will focus on light-element chemistry, with no
metals or other heavy atoms except at most one when scientifically essential,
because this is the relevant space for sustainable, energy-conscious molecular
design and because it keeps the physics interpretable.

The scientific challenge is therefore a scale jump: to extend toward 100 atoms
the kind of spectroscopic accuracy that today is realistically achieved only
for systems of about 10 atoms. MATRIX will not claim that the same protocol
should be pushed blindly to 1000 atoms. At that scale, and in condensed phases,
spectral resolution and molecular specificity usually decrease; the appropriate
description becomes multiscale, combining QM/MM, coarse-grained and continuum
models. MATRIX will supply the high-resolution molecular physics and
microsolvation-derived parameters needed by those larger-scale descriptions.

The second fragmentation is operational. Advanced spectroscopic workflows are
often locked into expert scripts, local conventions and fragile transfers
between quantum-chemistry outputs, visualization programs, databases and
analysis codes. MATRIX will make the same molecular state accessible through
command-line, graphical and web-aware interfaces, with execution strategies
tuned to the machine being used while preserving the same scientific model and
provenance.

This interoperability is planned from the beginning, not added after separate
tools have already diverged. Shared parsers, topology, symmetry, molecular
state, numerical kernels and validation layers allow each tool to be developed
scientifically one at a time while reusing the same assumptions and data.

## Central Hypothesis

A single enriched molecular state, carrying geometry, topology, symmetry,
fragment descriptors, internal coordinates, force fields, quantum-mechanical
data and provenance, can support a new class of predictive spectroscopy and
structure tools. If all observables are computed from this shared state, then
multi-level models can be combined without losing mode identity, symmetry
classification or chemical interpretability.

Systematic errors should first be removed by knowledge-based models: improved
coordinates, topology, symmetry, physically meaningful force constants,
electrostatics, dispersion, multi-level quantum data and different accuracy
levels for intra- and inter-fragment contributions. Only the remaining
residual, effectively random component should be learned. MATRIX will therefore
use machine learning as DeltaML, not as a black-box replacement for molecular
theory.

Accuracy will be treated by decomposing observables into contributions rather
than by applying one uniform level or neglecting expensive terms. Contributions
that are costly but small in the final observable may be computed at lower
level, provided that they are included coherently and can later be promoted to
higher level when needed. This generalizes the philosophy of
electronic-structure PCS/composite strategies and three-level vibrational
treatments to structures, force fields, spectra and fragments.

MATRIX is not only a technological or implementation advance. Some of its
components are new scientific methods: localized and symmetry-adapted
generalized internal coordinates, chemically controlled primitive-coordinate
selection, electronic-structure PCS/composite strategies, VPT2/VCI protocols
tied to persistent molecular state, and local-mode/fragment-aware models for
molecules beyond the current high-accuracy size regime. The ecosystem is
needed to make these methods usable, reproducible and comparable, but the
methodological advances are central to the ERC proposal.

## Breakthrough Objectives

1. Build a unified molecular state representation linking geometry, topology,
   symmetry, synthons, GICs and quantum observables.
2. Develop new symmetry- and fragment-aware internal coordinates for robust
   force fields, spectroscopy and refinement.
3. Develop original localized/symmetric GIC, electronic-structure
   PCS/composite, VPT2/VCI and local-mode methodologies for predictive
   spectroscopy beyond the current small-molecule accuracy frontier.
4. Generalize PCS and three-level vibrational logic into an accuracy hierarchy
   where all contributions are included coherently and can be upgraded
   independently.
5. Treat strongly correlated intra-fragment and weakly correlated
   inter-fragment contributions at different levels for energies, structures
   and force fields.
6. Create hybrid multi-level vibrational/rovibrational models that preserve
   normal-mode identity and chemical interpretation.
7. Introduce DeltaML corrections for residual errors only after systematic
   physical and chemical errors have been removed by the MATRIX model.
8. Extend the framework to frontier spectroscopies: vibronic spectra, VSCF
   local-mode anharmonicity, ionization potentials and electron affinities.
9. Define hand-off protocols from high-resolution MATRIX calculations to
   multiscale QM/MM, coarse-grained and continuum descriptions when molecular
   size or condensed-phase complexity makes full spectroscopic resolution
   unnecessary or impossible.
10. Create a robust, flexible and user-friendly ecosystem that communicates with
   external databases, web services, visualization tools, quantum-chemistry
   engines and remote computing resources.
11. Deliver reproducible prediction-to-publication workflows with traceable
   provenance and direct comparison to experimental databases.

## High-Risk / High-Gain Angle

The risk is that no single representation can remain chemically interpretable,
numerically stable and general enough across all relevant molecular regimes.
The gain is a step change in computational spectroscopy: workflows that can
combine quantum chemistry, internal-coordinate physics, symmetry, topology,
experimental databases and publication-grade analysis without uncontrolled
manual translation between tools.

The deliberately constrained chemical domain, light-element systems up to about
100 atoms, makes the project ambitious but not diffuse. The high-risk element
is to prove that knowledge-based error removal can extend near-small-molecule
accuracy to this frontier size range, and can reduce the residual error to a
regime where DeltaML is transferable and scientifically interpretable.
The gain is an accuracy model that is both economical and complete: no
contribution is silently discarded, but computational effort is assigned where
it has the largest effect on the final observable.

## Why Now

- Quantum-chemistry outputs are rich but underused because parsers and
  downstream tools are fragmented.
- Modern Python/Fortran hybrid workflows make it possible to combine robust
  legacy numerical kernels with reproducible high-level orchestration.
- Experimental spectra and molecular databases can now be integrated directly
  into automated validation loops.
- DeltaML is mature enough to be useful, but only if applied after physical
  sources of systematic error have been explicitly controlled.
- Hardware-aware workflows can now make the same scientific protocol usable on
  laptops, workstations and high-performance computing facilities.
- ERC-scale support can turn MATRIX from a strong research code into a
  scientific programme with new theory, algorithms and benchmarked applications.

## Provisional B1 Structure

1. The scientific bottleneck: molecular information is lost between tools.
2. The MATRIX idea: one molecular state, many rigorously coupled observables.
3. Methodological breakthroughs.
4. Representative scientific targets and demonstrators.
5. Why the PI and group can deliver this: prior record of coordinating theory,
   algorithms, implementation, validation and experiment in molecular
   spectroscopy.
