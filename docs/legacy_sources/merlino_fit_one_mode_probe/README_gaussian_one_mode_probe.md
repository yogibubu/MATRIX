# Gaussian One-Mode Probe

File:
- `/Users/vincenzobarone/merlino_fit/scripts/gaussian_one_mode_probe.py`

## Purpose

This script builds a **single-mode 1D rovibrational probe** from a Gaussian
anharmonic `log`/`fchk` pair.

It is designed for cases where one wants to inspect, for a selected normal
mode:

- the Gaussian-consistent normal coordinate `Q`,
- the corresponding inertia profile `I(Q)`,
- the inverse-inertia profile `1/I(Q)`,
- the 1D quartic potential built from Gaussian reduced force constants.

The script is intentionally conservative: it only performs the steps that are
already numerically validated.

## What is validated

For the selected mode, the script:

1. aligns the harmonic model to the Gaussian `dIdQ` convention;
2. reconstructs the mode displacement in the same coordinate convention;
3. verifies that finite-difference `dI/dQ` matches the `Inertia Moments Derivatives w.r.t. Normal Modes`
   block printed by Gaussian.

This means that the script provides a reliable Gaussian-consistent mapping:

- `Q -> geometry(Q)`
- `Q -> I(Q)`

## Potential convention

The script uses the 1D potential

\[
V(Q) = \frac12 \omega^2 Q^2 + \frac16 \phi_{iii} Q^3 + \frac1{24} \phi_{iiii} Q^4
\]

with:

- `omega` from Gaussian anharmonic force data,
- `phi_iii`, `phi_iiii` taken as **Gaussian reduced force constants**,
- `Q` in the same aligned coordinate convention used for the `dIdQ` validation.

The reduced-force-constant convention is the one already encoded in:

- `/Users/vincenzobarone/centrifugal/gaussian_force_constant_units.py`

## Important limitation

The script now also provides a **1D variational vs perturbative comparison**
for a selected rotational observable `B_tau(Q)`, using:

- a harmonic-oscillator basis,
- the quartic/cubic 1D Hamiltonian built from Gaussian reduced force constants,
- a polynomial representation of `B_tau(Q)`.

This is suitable for internal quantitative tests and for comparing:

- harmonic average,
- perturbative average,
- variational average.

However, this should still be regarded as a controlled numerical probe, not yet
as a fully finalized production framework.

So, at present, it is safe to use it for:

- `Q`-dependent geometry sampling,
- `I(Q)` sampling,
- `1/I(Q)` sampling,
- inspection of the 1D quartic potential,
- variational/perturbative comparison of `B_tau(Q)` for a chosen axis.

It should still **not yet** be treated as a final production tool unless the
1D model is independently benchmarked case by case.

## Usage

Example:

```bash
python /Users/vincenzobarone/merlino_fit/scripts/gaussian_one_mode_probe.py \
  --log /Users/vincenzobarone/merlino_fit/alpha_resonances_data/formic_rDSD.log \
  --fchk /Users/vincenzobarone/merlino_fit/alpha_resonances_data/formic.fchk \
  --mode 10
```

JSON output:

```bash
python /Users/vincenzobarone/merlino_fit/scripts/gaussian_one_mode_probe.py \
  --log /Users/vincenzobarone/merlino_fit/alpha_resonances_data/formic_rDSD.log \
  --fchk /Users/vincenzobarone/merlino_fit/alpha_resonances_data/formic.fchk \
  --mode 10 \
  --json
```

Optional grid controls:

- `--qmax`
- `--nq`

Optional 1D comparison controls:

- `--axis A|B|C`
- `--basis-size`
- `--poly-degree`

## Output

The script reports:

- alignment choice (`representation`, axis flips, mode signs),
- Gaussian `dIdQ`,
- finite-difference `dIdQ`,
- norm of the difference,
- `omega`, `phi_iii`, `phi_iiii`,
- sampled `I(Q)` and `1/I(Q)` values on the selected grid,
- when `--axis` is used:
  - harmonic average,
  - perturbative average,
  - variational average,
  - `variational - perturbative`.

## Current status

This is a validated probe for:

- the Gaussian-consistent coordinate,
- the Gaussian-consistent `I(Q)` profile,
- an internal 1D variational/perturbative comparison for `B_tau(Q)`.
