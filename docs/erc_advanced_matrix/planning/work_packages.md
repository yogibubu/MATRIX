# Work Packages

This is the first MATRIX ERC work-package map. It is intentionally scientific,
not code-module based.

## WP1 - Unified Molecular State And QM Adapters

Goal:

Create the shared molecular state that preserves geometry, topology, symmetry,
normal modes, force fields, electronic data and provenance across all tools.

Key tasks:

- Harden external adapters for Gaussian, Molpro, MRCC and future PySCF.
- Define final schemas for QM tensors, electronic states, orbitals and spectra.
- Add run manifests for reproducibility.
- Ensure parsers, topology, symmetry, numerical kernels and validators are
  shared from the start by all tools.
- Build online help/manual generation from tool contracts.

Outputs:

- versioned MATRIX molecular-state specification;
- shared interoperability and reuse contracts for all tools;
- adapter benchmark suite;
- reproducibility manifest model.

## WP2 - Symmetry, Topology And Internal Coordinates

Goal:

Develop NEO as the mathematical core for non-redundant, symmetry-adapted and
fragment-aware internal coordinates.

Key tasks:

- Complete special coordinate theory for fragments, rings and virtual centres.
- Develop localized and symmetry-adapted GICs as new scientific coordinates,
  not only as Gaussian-compatible input constructs.
- Define chemically controlled primitive-coordinate selection strategies that
  preserve essential special coordinates.
- Maintain Python/Fortran parity.
- Connect topology/synthon descriptors to force-field and spectroscopic models.
- Define intra-fragment and inter-fragment level partitioning for structures
  and force fields.
- Prepare next-release nano-LEGO and fragment assembly foundations.

Outputs:

- NEO coordinate theory papers;
- primitive-coordinate selection and localized/symmetric GIC methodology
  papers;
- fragment-level force-field and structure partitioning protocols;
- benchmark corpus for GIC/SALC/B-matrix validation;
- fragment-aware coordinate protocols.

## WP3 - Predictive Vibrational And Rovibrational Spectroscopy

Goal:

Build rigorous multi-level vibrational and rotational spectroscopy workflows.

Key tasks:

- Extend hybrid level1+level2 models with mode-overlap validation.
- Generalize three-level vibrational logic so harmonic, anharmonic,
  equilibrium and vibrational-correction contributions can be computed at
  different coherent levels.
- Develop VSCF harmonically coupled anharmonic oscillators in local modes.
- Develop persistent-state VPT2/VCI protocols and integrate DVR in common
  analysis pipelines.
- Homologate WMS-Rot and connect rotational predictions to databases.

Outputs:

- predictive IR/Raman/VCD/ROA workflows;
- generalized three-level vibrational accuracy protocols;
- local-mode VSCF algorithms;
- persistent-state VPT2/VCI methods;
- rotational line-list validation protocols.

## WP4 - Electronic, Vibronic, IP And EA Workflows

Goal:

Extend MATRIX beyond ground-state vibrational spectroscopy to electronic and
charge-changing observables.

Key tasks:

- Add vibronic spectroscopy data contracts and workflows.
- Add ionization-potential and electron-affinity protocols.
- Add PySCF launcher and analysis support.
- Connect Molden/Avogadro/MOrbVis orbital visualization to reproducible state.

Outputs:

- vibronic spectra workflows;
- vertical/adiabatic IP/EA protocols;
- integrated orbital and electronic-state analysis.

## WP5 - Knowledge-Based Error Removal, DeltaML And Validation

Goal:

Make MATRIX predictions quantitatively reliable by separating systematic errors
that must be removed by physical and chemical knowledge from residual errors
that can be treated by DeltaML, and by transferring high-resolution molecular
information to multiscale models when size or phase makes full spectroscopic
resolution inappropriate.

Key tasks:

- Define the error hierarchy: coordinates, topology/symmetry, force constants,
  electrostatics, dispersion, anharmonicity, electronic-state model and
  residual statistical error.
- Generalize electronic-structure PCS/composite strategies to a
  contribution-based accuracy hierarchy across structures, force fields and
  spectra.
- Include computationally expensive but small contributions coherently at
  lower level instead of neglecting them, with explicit promotion paths to
  higher accuracy.
- Separate strongly correlated intra-fragment contributions from weaker,
  more transferable inter-fragment contributions.
- Build DeltaML models only on residuals after the knowledge-based correction
  layers have been applied.
- Restrict training and validation to chemically coherent light-element domains
  up to about 100 atoms.
- Define the scale hand-off from high-resolution MATRIX calculations to QM/MM,
  coarse-grained and continuum descriptions for 1000-atom systems and condensed
  phases.
- Derive transferable MM parameters from accurate microsolvation calculations.
- Add experimental spectrum database adapters.
- Build publication export backends for spectra and tables.
- Create curated benchmark corpora.
- Develop GUI workflows and in-line help for major tools.

Outputs:

- knowledge-based error decomposition protocols;
- contribution-based PCS/three-level accuracy protocols;
- intra-/inter-fragment accuracy partitioning models;
- DeltaML residual-correction models with uncertainty diagnostics;
- microsolvation-to-MM parameterization protocols;
- scale-bridging interfaces to QM/MM, coarse-grained and continuum models;
- benchmarked prediction-to-publication pipeline;
- database comparison tools;
- release-quality manuals and training material.

## WP6 - Robust Ecosystem And Adaptive Execution

Goal:

Turn MATRIX into a usable scientific ecosystem in which the same molecular
state and the same assumptions can move across command-line tools, graphical
interfaces, web/database resources and computing environments from laptops to
workstations and high-performance computing centres.

Key tasks:

- Define hardware-aware execution profiles for laptops, workstations and HPC
  centres.
- Add robust external interfaces to web services, experimental databases,
  molecular repositories, visualization tools and remote job launchers.
- Provide GUI workflows for spectroscopy, structure/synthons,
  thermo/kinetics, QM jobs, visualization and publication export.
- Ensure that every tool remains independently usable when the shared molecular
  state contains the required section.
- Preserve provenance, environment, software versions and execution manifests
  across local and remote runs.

Outputs:

- adaptive execution profiles and benchmark reports;
- web/database interface contracts;
- GUI workflows for major scientific tasks;
- reproducible local/workstation/HPC run manifests;
- user-facing manuals and online help for the main tools.
