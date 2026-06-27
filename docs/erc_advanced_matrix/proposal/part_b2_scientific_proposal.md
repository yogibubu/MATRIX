# Part B2 - Scientific Proposal Draft

Target length: follow the official ERC-2026-AdG template.

## Provisional Scientific Structure

1. State of the art and its limits.
2. MATRIX scientific concept.
3. Work Package 1 - Unified molecular state and adapters.
4. Work Package 2 - Symmetry/topology/GIC theory.
5. Work Package 3 - Predictive vibrational and rovibrational spectroscopy.
6. Work Package 4 - Electronic and vibronic observables.
7. Work Package 5 - Knowledge-based error removal, DeltaML, validation and
   publication-grade workflows.
8. Work Package 6 - Robust ecosystem, web interfaces and adaptive execution.
9. Risk management and alternative strategies.
10. Resources, team and lump-sum justification.

## Core Technical Claims To Develop

- The enriched `xyzin` state is not just a file format; it is a scientific
  representation preserving molecular identity across calculations.
- Interoperability and tool reuse are design principles from the beginning, so
  each scientific module can advance independently while sharing parsers,
  topology, symmetry, numerical kernels and validation.
- NEO/GICs provide chemically interpretable coordinates for force fields,
  spectra, refinement and optimization.
- Localized and symmetry-adapted GICs, electronic-structure PCS/composite
  strategies, persistent-state VPT2/VCI and local-mode models are original
  scientific methods, not only implementation features.
- Accuracy is decomposed contribution by contribution: expensive terms with
  limited final weight are computed at coherent lower levels rather than
  neglected, and each contribution remains independently upgradeable.
- Intra-fragment and inter-fragment effects can use different levels not only
  for energies but also for structures and force fields.
- Hybrid level1+level2 spectra become scientifically safe only when normal
  modes are explicitly matched.
- Fragment/topology/synthon descriptors allow chemically meaningful transfer of
  parameters, scaling and uncertainty.
- Machine learning is used only as DeltaML on residual random-like errors after
  systematic errors have been removed by knowledge-based molecular physics.
- Direct links to experimental databases convert MATRIX into an inference and
  validation framework, not only a calculation runner.
- Hardware-aware execution makes one scientific protocol usable on laptops,
  workstations and high-performance computing centres.
- Web/database interfaces connect MATRIX to experimental spectra, molecular
  repositories, visualization tools, remote job launchers and external
  scientific services.

## Target Domain

MATRIX will target molecular systems up to about 100 atoms. The goal is to
bring to this range the kind of spectroscopic accuracy that is currently
reachable only for much smaller systems, roughly up to 10 atoms. This range
defines the current frontier where spectroscopy, quantum chemistry and
molecular nano-engineering meet: systems are large enough for real function,
conformational richness and weak interactions, but still accessible to
interpretable high-resolution modelling.

The proposal should focus on metal-free and heavy-atom-light chemistry, allowing
at most one heavy atom when scientifically essential. This boundary strengthens
the scientific message: MATRIX is designed for sustainable, abundant-element
molecular design and avoids making relativistic or heavy-element effects the
main source of complexity.

The proposal should make the scale logic explicit. MATRIX does not aim to treat
1000-atom systems with the same high-resolution isolated-molecule protocol,
because at that size spectral resolution and molecular specificity normally
decrease. For those systems, and for condensed phases, the correct level is a
multiscale one: QM/MM, coarse-graining and continuum models. MATRIX should
generate the accurate molecular building blocks for those models, including MM
parameters from microsolvated clusters treated at high accuracy.

## Demonstrator Classes

Choose 3-5 classes that show frontier value:

- polycyclic and fused-ring systems for GIC/ring-puckering complexity;
- weakly bound complexes and hydrogen-bonded systems for fragment-centre
  coordinates and low-frequency modes;
- non-covalent systems, ring complexes and single-heavy-atom edge cases where
  atom-to-ring or atom-to-bond centres matter;
- fluxional systems for DVR/VSCF/local-mode anharmonicity;
- molecules with well-characterised rotational/vibrational/electronic spectra.
- light-element molecular nano-engineering targets approaching 100 atoms,
  where spectroscopy remains measurable but current modelling loses robustness.

## Original Methodological Advances

The proposal must make clear that MATRIX is not justified only by better
software engineering. The ERC-level novelty includes:

- localized and symmetry-adapted generalized internal coordinates with
  chemically interpretable reduction and analytic B matrices;
- chemically controlled primitive-coordinate selection that preserves
  essential coordinates while removing uncontrolled redundancy;
- electronic-structure PCS/composite strategies generalized into a broader
  accuracy hierarchy for observables;
- VPT2/VCI workflows linked to a persistent molecular state, so anharmonic
  treatments inherit topology, symmetry, mode identity and provenance;
- local-mode and fragment-aware models for systems larger than the present
  high-accuracy regime;
- controlled multi-level corrections in which mode identity and force-field
  meaning are preserved before any DeltaML residual correction is attempted.
- intra-/inter-fragment level partitioning for energies, equilibrium
  structures, vibrational corrections and force fields.

## Expected Outcomes

- New theory and algorithms for reusable molecular state propagation.
- New internal-coordinate and symmetry-adaptation methods.
- New electronic-structure PCS/composite, VPT2/VCI and local-mode
  methodologies embedded in a unified molecular-state framework.
- Generalized contribution-based accuracy protocols for electronic,
  structural, vibrational and fragment observables.
- New multi-level vibrational and rovibrational prediction protocols.
- A controlled DeltaML layer for residual errors after knowledge-based
  correction of systematic effects.
- Scale-bridging protocols to transfer high-resolution molecular information
  into QM/MM, coarse-grained and continuum descriptions.
- A robust command-line, GUI and web-aware ecosystem tuned to the available
  computing environment.
- New tools for vibronic, IP/EA and local-mode anharmonic spectroscopy.
- Curated benchmark corpora and reproducible workflows.
- Publication-ready spectra/tables linked to full computational provenance.
