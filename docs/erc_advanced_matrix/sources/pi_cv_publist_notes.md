# PI CV And Publication-List Notes

Sources inspected:

- `/Users/vincenzobarone/Desktop/CV_VB copia.pdf`
- `/Users/vincenzobarone/Desktop/publist_VB/publist.tex`
- `/Users/vincenzobarone/Desktop/publist_VB/VB.bib`
- `/Users/vincenzobarone/Desktop/publist_VB/verified_200_citations.csv`
- `/Users/vincenzobarone/Desktop/publist_VB/pubs_per_year_chart.tex`
- `/Users/vincenzobarone/Desktop/publist_VB/pubs_2021_2025_chart.tex`
- `/Users/vincenzobarone/Desktop/publist_VB/pubs_per_decade_pie.tex`
- `/Users/vincenzobarone/Desktop/publist_VB/citations_rank_chart.tex`
- corresponding PDF graphics in `/Users/vincenzobarone/Desktop/publist_VB/`

These notes are working evidence for the ERC B1 CV/track-record section. They
should be checked against the final official ERC CV template before submission.

## CV Facts To Use

- International profile in theoretical and computational chemistry.
- Major contributions named in the CV: PBE0, CPCM/continuum solvation models,
  vibrational and vibronic spectroscopy, Gaussian development and Pisa
  Composite Schemes.
- Handbook of Electronic Structure Theory, Elsevier, 2026, co-edited with
  Majdi Hochlaf, over 700 pages.
- Continued scientific activity after retirement in 2023, with the CV stating
  more than 60 publications since then.
- About 950 peer-reviewed journal publications, about 96000 citations, more
  than 88000 citations excluding self-citations, H-index 105.
- CV states 11 publications above 1000 citations; the current top-citation CSV
  contains 12 entries at or above 1000 citations and should be reconciled with
  the final citation source/date.
- More than 150 invited lectures and more than 100 plenary lectures.
- Full professor at Scuola Normale Superiore, 2009-2023; full professor at
  Universita di Napoli Federico II, 1994-2008.
- PI of national and international projects, including ERC Advanced Grant
  DREAMS, over EUR 2M.
- Founder/promoter of SMART Lab and the AVOGADRO data center.
- Coordinator of the INSTM Skies Village reference center, linking groups in
  Milan, Bologna and Pisa.
- Accademia dei Lincei, Academia Europaea, Fellow of the Royal Society of
  Chemistry, International Academy of Quantum Molecular Sciences.
- Director/Rector of Scuola Normale Superiore, 2016-2019.
- President of Societa Chimica Italiana; director of CNR research institutes.
- Member of CNGR and President of the Chemistry panel for VQR/ANVUR; evaluator
  for agencies including NSF, ANR, FWF and FNRS.

## Publication-List Facts

The `publist_VB` directory already contains statistics and graphics. Use these
files as the primary internal source for production/citation plots rather than
recomputing them in the proposal:

- `pubs_per_year_chart.tex/pdf`: yearly publication count from 1977 to 2026.
- `pubs_2021_2025_chart.tex/pdf`: recent yearly publication count, 2021-2026.
- `pubs_per_decade_pie.tex/pdf`: decade-grouped publication distribution.
- `citations_rank_chart.tex/pdf`: citation-bin distribution for the top 200
  verified works.
- `verified_200_citations.tex/pdf/csv`: tabular source for the top 200
  citation records.
- `publist.pdf`: full publication list with embedded graphics.

The existing decade chart reports `Total considered: 923`.

The existing 2021-2026 chart reports:

  - 2021: 24
  - 2022: 22
  - 2023: 23
  - 2024: 18
  - 2025: 24
  - 2026: 16

The existing decade chart reports:

- 1976-1985: 72 publications, 7.8%
- 1986-1995: 99 publications, 10.7%
- 1996-2005: 178 publications, 19.3%
- 2006-2015: 317 publications, 34.3%
- 2016-2025: 241 publications, 26.1%
- 2026-2035: 16 publications, 1.7%

The existing citation-bin chart for the top 200 verified works reports:

- 50-99 citations: 88 publications
- 100-199 citations: 64 publications
- 200-499 citations: 30 publications
- 500-999 citations: 6 publications
- 1000-1999 citations: 4 publications
- 2000-4999 citations: 5 publications
- 5000+ citations: 3 publications

Independent consistency checks performed during this note preparation:

- `publist.tex` contains 923 `pubitem` entries across 1977-2026, matching the
  decade chart total.
- `publist.tex` contains 127 entries for 2021-2026.
- `publist.tex` contains 81 entries for 2023-2026.
- `verified_200_citations.csv` lists 200 highly cited works; among valid rows,
  12 are at or above 1000 citations, 18 at or above 500 citations, and 48 at or
  above 200 citations. This should be reconciled with the CV's current statement
  of 11 publications above 1000 citations, probably by fixing the source/date
  used for the final CV.

## MATRIX-Relevant Publication Evidence

Use these as examples, not as the final selected publication list.

| MATRIX claim | Evidence from publication list | Relevance |
| --- | --- | --- |
| Foundational electronic-structure methods | "Toward Reliable Density Functional Methods without Adjustable Parameters: The PBE0 Model", J. Chem. Phys., 1999, DOI 10.1063/1.478522 | Demonstrates ability to introduce methods adopted worldwide. |
| Solvation and environment | "Quantum calculation of molecular energies and energy gradients in solution by a conductor solvent model", J. Phys. Chem. A, 1998, DOI 10.1021/jp9716997; CPCM work in J. Comput. Chem., 2003, DOI 10.1002/jcc.10189 | Supports MATRIX treatment of environment, microsolvation and condensed-phase hand-off. |
| Automated anharmonic spectroscopy | "Anharmonic vibrational properties by a fully automated second order perturbative approach", J. Chem. Phys., 2005, DOI 10.1063/1.1824881 | Directly de-risks VPT2/GVPT2 and large-scale anharmonic workflows. |
| Vibronic spectroscopy | "Fully integrated approach to compute vibrationally resolved optical spectra: from small molecules to macrosystems", J. Chem. Theory Comput., 2009, DOI 10.1021/ct8004744; "General time dependent approach to vibronic spectroscopy including Franck-Condon, Herzberg-Teller, and Dushinsky effects", J. Chem. Theory Comput., 2013, DOI 10.1021/ct400450k | Supports vibronic and electronic-spectroscopy objectives. |
| Spectroscopy review leadership | "Accuracy and Interpretability: The Devil and the Holy Grail. New Routes across Old Boundaries in Computational Spectroscopy", Chem. Rev., 2019, DOI 10.1021/acs.chemrev.9b00007; "Computational Molecular Spectroscopy", Nat. Rev. Methods Primers, 2021, DOI 10.1038/s43586-021-00034-1 | Establishes the conceptual frame of accuracy plus interpretability. |
| Theory-experiment structural work | "Semi-experimental equilibrium structure determinations by employing B3LYP/SNSD anharmonic force fields", J. Phys. Chem. A, 2015, DOI 10.1021/jp511432m | Supports MORPHEUS/SEFit and the theory-experiment feedback loop. |
| PCS and accuracy hierarchy | PCS papers from 2023-2026, including "PCS/Bonds", J. Chem. Phys., 2023, DOI 10.1063/5.0167296; "Accurate Thermochemical and Kinetic Parameters at Affordable Cost by Means of the Pisa Composite Scheme", J. Chem. Theory Comput., 2023, DOI 10.1021/acs.jctc.3c00817; "From Energies to Structures: PNO-F12 Geometry Optimization of Large Molecules for Near-Spectroscopic Accuracy", J. Chem. Theory Comput., 2026, DOI 10.1021/acs.jctc.6c00317 | Direct precedent for the contribution-based accuracy strategy. |
| Synthons, fragments and NanoLEGO | "Toward Spectroscopic Accuracy for the Structures of Large Molecules at DFT Cost: Refinement and Extension of the Nano-LEGO Approach", J. Phys. Chem. A, 2023, DOI 10.1021/acs.jpca.3c01617; "Accurate geometries of Large Molecules by integration of the Pisa Composite Scheme and the Templating Synthon Approach", J. Phys. Chem. A, 2024, DOI 10.1021/acs.jpca.3c08382; "High-Fidelity Ring Fragments for Molecular Design and Spectroscopy: the PCS-LCB25-Nano-LEGO Framework", J. Chem. Theory Comput., 2025, DOI 10.1021/acs.jctc.5c01188 | De-risks the fragment/synthon/nano-LEGO future release and supports the 100-atom frontier. |
| Rotational spectroscopy and refinement | "From Quantum-Chemical Predictions to Rotational Spectral Assignment and Refinement", J. Chem. Phys., 2026, DOI 10.1063/5.0339084 | Supports the rotational-spectroscopy branch and theory-experiment integration. |
| Anharmonicity in extended environments | "Anharmonicity in Molecular Crystals: Generalized Perturbation Theory Meets Periodic Computations", J. Phys. Chem. Lett., 2025, DOI 10.1021/acs.jpclett.5c02217 | Supports the hand-off between high-resolution molecular models and condensed-phase/periodic contexts. |

## ERC Argument

The CV and publication list support a strong ERC Advanced Grant narrative:

- The PI has already produced globally adopted theory and computational tools.
- The PI has repeatedly connected theory, implementation and experiment.
- The PI has coordinated complex scientific and institutional structures.
- The PI's post-2023 activity shows current scientific momentum and can be
  framed as increased scientific focus after removal of administrative duties.
- Recent post-retirement work includes ideas directly aligned with MATRIX:
  PCS extensions, NanoLEGO/synthons/LCB25, accuracy ladders,
  WMSRot/refinement and local-correlation-enhanced structures.
- MATRIX is a natural next step, not a late-career change of direction.
