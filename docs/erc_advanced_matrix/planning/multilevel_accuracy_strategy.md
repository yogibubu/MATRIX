# Multi-Level Accuracy Strategy

## Principle

MATRIX should pursue accuracy by controlled decomposition, not by using one
uniform expensive method for all terms and not by neglecting small terms. Every
physically meaningful contribution should remain in the model. Contributions
that are computationally expensive but have limited weight in the final
observable may be computed at lower level, provided that the approximation is
coherent, traceable and independently upgradeable.

## Relation To Existing Ideas

This generalizes two existing philosophies:

- electronic-structure PCS/composite strategies, where different energy
  contributions can be evaluated at different levels;
- three-level vibrational strategies, where harmonic, anharmonic, equilibrium
  and vibrational-correction terms can be assigned different accuracies.

MATRIX should make this logic general across electronic structure, equilibrium
geometries, vibrational corrections, force fields, spectra and fragment
interactions.

## Fragment-Level Accuracy

Intra-fragment and inter-fragment interactions should be treated differently.
Intra-fragment terms are often strongly correlated and chemically specific;
inter-fragment terms are usually weaker, more transferable and more local.
Existing local-correlation methods exploit this distinction mainly for
energies. MATRIX should extend the same idea to:

- equilibrium structures;
- vibrational corrections;
- force constants and force fields;
- fragment-centre coordinates and low-frequency modes;
- non-covalent and long-range interactions.

## Upgrade Path

The molecular state should record the level used for each contribution. If a
contribution becomes important for a given observable, it can be recomputed at
higher level without changing the parser, topology, symmetry, coordinate
definitions or downstream analysis.

## DeltaML Layer

DeltaML is the final layer. It should learn only the residual random-like error
that remains after knowledge-based decomposition, multi-level treatment,
fragment partitioning and physical corrections have been applied.

## ERC Message

The novelty is not simply that MATRIX automates multi-level calculations. The
scientific point is a general accuracy calculus for spectroscopy and molecular
properties: include all contributions, assign effort where it matters, preserve
provenance and make every approximation reversible.
