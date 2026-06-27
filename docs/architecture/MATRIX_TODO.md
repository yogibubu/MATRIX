# MATRIX TODO

## WMS-Rot Homologation

- Remove the hard `pandas` dependency from MATRIX-facing WMS-Rot services after
  the line-list, assignment and CSV layers have typed ORACLE records and golden
  regression tests.
- Keep `pandas` available only as a temporary compatibility dependency for the
  imported first-party WMS-Rot engine.
- Move WMS-Rot adapter boundaries toward ordinary MATRIX services: shared
  xyzin sections, run manifests, publication exports, shared diagonalizer and
  shared GUI workflow state.
- Do not add Gaussian, QM-output or topology parsing inside WMS-Rot. All such
  data must enter through shared MATRIX adapters and enriched xyzin sections.

## DVR Diagonalization

- Keep large DVR Hamiltonian diagonalizations behind
  `oracle_core.diagonalizer`.
- Do not add private DVR diagonalizer wrappers in new Python DVR workflows.
- Add GPU/performance regression cases when realistic large DVR fixtures are
  available.

## Experimental Spectrum Databases

- Add a shared experimental-spectrum adapter layer for external databases
  instead of hardcoding database queries inside spectroscopy GUIs.
- Vibrational: NIST IR data may be fetched automatically only for gas-phase
  JCAMP records; condensed-phase and missing records require user instruction.
- Rotational: add database-backed comparison for predicted line lists against
  gas-phase microwave/submillimeter catalogs such as CDMS/JPL/Splatalogue when
  a reliable molecule identifier and frequency/unit contract are available.
- Keep database provenance, query URL, phase/state and selection filters in the
  export metadata so publication plots remain traceable.
