# MATRIX ERC Advanced Grant 2026

This directory is the working container for the ERC Advanced Grant proposal
built around MATRIX as the scientific project.

Working title:

> MATRIX: a unified molecular-intelligence framework for predictive
> spectroscopy, structure and dynamics

The title is provisional. The final ERC title should be shorter, more
frontier-research oriented and less software-centred.

## Current Call Facts To Recheck

- Call: ERC-2026-AdG.
- Official deadline shown on the ERC Advanced Grant page: 27 August 2026.
- Scheme: Advanced Grant for established, leading PIs.
- Funding: up to EUR 2.5M for up to 5 years, plus eligible additional funding.
- ERC criterion: excellence as the sole criterion.
- 2026 Advanced Grants are implemented as a lump-sum pilot.

Always re-check the official ERC page, Work Programme 2026, Information for
Applicants 2026 and Funding & Tenders Portal before freezing the proposal.

## Directory Map

| Directory | Purpose |
| --- | --- |
| `proposal/` | Scientific text drafts: B1, CV/track record, B2. |
| `planning/` | Work packages, timeline, budget logic, risks, panel strategy. |
| `admin/` | Host institution, ethics, open science and submission checklist. |
| `figures/` | Figure concepts and later source files for proposal figures. |
| `reviews/` | Internal review log, panel-reader comments and revision history. |
| `sources/` | Official ERC references and call notes. |

## Immediate Writing Priorities

1. Define the ERC-level scientific question in one paragraph, not as a software
   refactoring story.
2. Identify the target ERC panel, likely within Physical Sciences and
   Engineering, and map keywords to panel descriptors.
3. Write a one-page high-risk/high-gain concept note before expanding B1.
4. Convert MATRIX components into scientific objectives and hypotheses.
5. Build a credible five-year work plan with milestones, risks and fallback
   strategies.
6. Start budget design early because ERC-2026-AdG uses a lump-sum model.

## First Boundary

For the ERC proposal, MATRIX should not be presented as "a code base". It
should be presented as a scientific platform enabling new chemistry:

- unified molecular representations across quantum chemistry, spectroscopy,
  thermochemistry and dynamics;
- rigorous propagation of geometry, topology, symmetry and uncertainty through
  all molecular observables;
- hybrid multi-level vibrational and rovibrational models;
- predictive spectra and molecular descriptors with reproducible provenance;
- extensible frontier modules for vibronic spectroscopy, VSCF local modes,
  PySCF-driven electronic structure, ionization potentials and electron
  affinities.

## Scientific Scope

MATRIX should be positioned on molecular systems up to about 100 atoms, aiming
to bring to this size range the spectroscopic accuracy that is currently
reachable only for systems of roughly 10 atoms. This is large enough to be a
current frontier for high-resolution spectroscopy, predictive nano-engineering
and realistic molecular function, but still small enough for physically
interpretable quantum-chemical and spectroscopic models.

The target chemical domain is light-element molecular science: organic,
bio-inspired, atmospheric, supramolecular and nano-engineered systems without
metals or other heavy atoms, with at most one heavy atom when scientifically
essential. This boundary is scientifically useful and strategically coherent
with energy and environmental arguments: MATRIX should emphasize sustainable,
abundant-element molecular design rather than heavy-element chemistry.

Machine learning is part of the project only as DeltaML for residual,
effectively random errors after systematic errors have been removed by
knowledge-based physics, chemistry, topology, symmetry and internal-coordinate
models. This distinction is essential for the ERC narrative: MATRIX does not
replace chemical understanding with black-box learning; it uses learning only
where the remaining error is no longer systematically explainable by the
knowledge-based model.

The proposal should explicitly explain why the target is not 1000 atoms. At
that scale, and in condensed phases, spectral resolution and molecular
specificity generally decrease, and the correct description becomes multiscale:
QM/MM, coarse-graining and continuum models. MATRIX should provide the accurate
small- and medium-size molecular physics that feeds those models, including MM
parameters derived from accurate microsolvation calculations, rather than
pretending that a high-resolution isolated-molecule protocol remains the right
tool at every scale.

## Ecosystem Boundary

MATRIX should also be presented as a robust scientific ecosystem. The
knowledge-based theory is the core, but the ERC project must make it usable
across real research settings: laptops for inspection and teaching,
workstations for routine production, and high-performance computing centres
for large benchmark and production campaigns.

The ecosystem should be flexible and user-friendly without hiding the physics.
It should communicate with external resources through documented web/database
interfaces, including experimental spectra, molecular databases, web services,
remote job launchers and visualization tools. The central design principle is
that the same molecular state and the same scientific assumptions travel across
all environments, while execution is tuned to the available machine.

Interoperability is not an add-on. It is part of the scientific design from the
start: each tool should reuse the same parsers, topology, symmetry, molecular
state, numerical kernels and validation layers. This makes it possible to
develop one scientific tool at a time without rebuilding infrastructure or
changing assumptions between modules.

## Accuracy Philosophy

MATRIX should treat accuracy as a controlled decomposition of contributions.
No physically meaningful contribution should be ignored, but computationally
expensive terms that contribute only a limited percentage to the final
observable may be computed at a lower level, provided that this is done
coherently, traceably and with an explicit upgrade path. When higher accuracy
is needed, any contribution can be promoted to a higher level without changing
the molecular state or the rest of the workflow.

This generalizes the philosophy behind electronic-structure PCS/composite
strategies and three-level vibrational treatments: equilibrium structures and
vibrational corrections, harmonic and anharmonic terms, local and non-local
contributions can be computed at different levels, but all are kept in the
model. Strongly correlated intra-fragment interactions and more weakly
correlated inter-fragment interactions should be allowed to use different
levels not only for energies, as in local-correlation methods, but also for
structures and force fields. DeltaML is then the final layer for residual
random-like errors after this knowledge-based decomposition has been applied.
