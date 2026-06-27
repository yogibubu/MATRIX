# ORACLE GUI Architecture

The ORACLE GUI is a project dashboard over the shared enriched XYZ container.
It must not implement private chemistry, parser, fitting or Fortran logic.

The allowed call path is:

```text
GUI view -> oracle_gui controller -> oracle-* service/CLI -> xyzin sections
```

The forbidden call paths are:

```text
GUI view -> private XYZ/Gaussian/Z-matrix parser
GUI view -> private topology or synthon perception
GUI view -> direct Fortran executable without an oracle-* service boundary
```

## Windows

The first window is `ORACLE Project Dashboard`. It opens or creates the active
project, displays the molecule state and lists all known `xyzin` sections.
Every other window is reachable from this state.

The planned windows are:

| Window | Responsibility | Main sections |
| --- | --- | --- |
| ORACLE Project Dashboard | Project state, validation, section/workflow status | all sections |
| ORACLE-Babel / Preprocessing | Import XYZ, QM formats, Z-matrix, SMILES/RDKit and LCB25 | `#SOURCE`, `#BASIC`, `#SYMMETRY`, `#TOPOLOGY`, `#SYNTHONS` |
| Molecule Editor / Avogadro Bridge | Open/edit the first XYZ block and reimport edited coordinates | XYZ block |
| Topology / Synthons / Fragments | Inspect bonds, rings, charges, synthons and fragments | `#TOPOLOGY`, `#SYNTHONS`, `#FRAGMENTS` |
| GICForge | Build, symmetrize and diagnose GICs and B matrices | `#GIC`, `#SYCART` |
| GF / PED | Harmonic force-field analysis from Hessian plus GICs | `#CARTESIAN_HESSIAN`, `#GF_PED` |
| SEFit / MORPHEUS | Single-molecule and ensemble semiexperimental refinement | `#ISOTOPOLOGUES`, `#MORPHEUS` |
| Thermo / Rovib | Rotational/vibrational summaries and thermochemistry | `#ROTATIONAL`, `#VIBRATIONAL`, `#THERMO` |
| Anharmonic: VPT2 / VCI / DVR | Run and collect anharmonic workflow state | `#QFF`, `#VPT2_VCI`, `#DVR` |
| QM Jobs | Generate Gaussian inputs and normalize QM output sections | `#CARTESIAN_HESSIAN`, `#NORMAL_MODES`, `#QFF` |
| Diagnostics / Regression | Corpus, Python/Fortran and benchmark audits | reports/artifacts |

## Implementation Boundary

`oracle-gui` owns only:

- project and section view-models;
- GUI command specifications;
- optional Qt widgets;
- file selection and user interaction state.

Scientific behavior remains in the owning packages:

- `oracle-chem` for ORACLE-Babel, symmetry, topology and synthons;
- `oracle-gicforge` for GIC construction, symmetrization and B matrices;
- `oracle-gf` for GF/PED;
- `oracle-morpheus` for SEFit and MORPHEUS;
- `oracle-rovib`, `oracle-thermo`, `oracle-vpt2-vci` and `oracle-dvr` for
  their corresponding sections.

The source-tree entry point is:

```bash
python -m oracle_gui [molecule.xyzin]
```

The installable console entry point is:

```bash
oracle-gui [molecule.xyzin]
```

Qt is optional and loaded only when the GUI process starts. Headless tests cover
the controllers without requiring a display server.
