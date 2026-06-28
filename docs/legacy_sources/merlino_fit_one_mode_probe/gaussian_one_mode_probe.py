#!/usr/bin/env python3
"""One-mode Gaussian probe for rovibrational diagnostics.

This utility does the parts that are now under control for a single Gaussian
anharmonic log/fchk pair:

1. align the harmonic model to Gaussian's printed dIdQ convention;
2. extract one selected normal mode in that convention;
3. build I(Q) and 1/I(Q) by Cartesian displacement along the aligned mode;
4. build the 1D quartic potential using Gaussian reduced force constants.

The script is intentionally conservative: it reports and samples the 1D model,
but it does not yet claim a fully validated 1D Schrödinger solver in the same
coordinate convention. The output is therefore suitable for analysis,
inspection, and later reuse in a validated solver.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


CENTRIFUGAL = Path("/Users/vincenzobarone/centrifugal")
if str(CENTRIFUGAL) not in sys.path:
    sys.path.insert(0, str(CENTRIFUGAL))

from gaussian_vpt_parser import (  # noqa: E402
    parse_gaussian_anharmonic_force_data,
    parse_gaussian_fchk_harmonic_data,
    parse_gaussian_harmonic_data,
)
from rovib_distortion import harmonic_inertia_model_from_geometry_hessian  # noqa: E402
from vibrot_alpha import alpha_matrix_from_harmonic_and_cubic_cm  # noqa: E402


IDX_MAP = {(0, 0): 0, (1, 0): 1, (1, 1): 2, (2, 0): 3, (2, 1): 4, (2, 2): 5}
PAIR_LABELS = ("Ixx", "Ixy", "Iyy", "Ixz", "Iyz", "Izz")
AXES = ("A", "B", "C")
B_MHZ_FROM_I_AMU_A2 = 6.62607015e-34 / (8.0 * math.pi**2 * 1.66053906660e-27 * (1.0e-10) ** 2) / 1.0e6


@dataclass
class AlignmentResult:
    representation: str
    axis_flip: tuple[int, int, int]
    printed_to_harmonic: list[int]
    harmonic_to_printed: list[int]
    mode_signs: list[int]
    didq_error: float


def _parse_didq_block(log_path: str | Path) -> np.ndarray:
    lines = Path(log_path).read_text(encoding="utf-8").splitlines()
    start = [i for i, line in enumerate(lines) if "Inertia Moments Derivatives w.r.t. Normal Modes" in line][-1]
    rows = []
    for line in lines[start + 6 : start + 6 + 1000]:
        m = re.match(r"\s*Q\(\s*(\d+)\)\s+(.+)$", line)
        if m is None:
            break
        rows.append([float(x.replace("D", "E")) for x in m.group(2).split()[:6]])
    if not rows:
        raise ValueError("Could not parse Gaussian dIdQ block.")
    return np.array(rows, dtype=float).T


def _didq_from_model(model) -> np.ndarray:
    n_atoms = model.masses_amu.size
    modes = model.vib_vecs_mw_pa.reshape(n_atoms, 3, -1)
    didq = np.zeros((6, modes.shape[2]), dtype=float)
    ijx = 0
    for ix in range(3):
        for jx in range(ix + 1):
            for mode in range(modes.shape[2]):
                acc = 0.0
                for atom in range(n_atoms):
                    mass = math.sqrt(model.masses_amu[atom])
                    acc -= mass * model.coords_pa_ang[atom, ix] * modes[atom, jx, mode]
                    if ix == jx:
                        for kx in range(3):
                            acc += mass * model.coords_pa_ang[atom, kx] * modes[atom, kx, mode]
                didq[ijx, mode] = 2.0 * acc
            ijx += 1
    return didq


def _apply_mode_reorder_and_signs(model, order: np.ndarray, signs: np.ndarray) -> None:
    model.vib_vecs_mw_pa[:] = model.vib_vecs_mw_pa[:, order] * signs
    model.vib_freq_cm[:] = model.vib_freq_cm[order]


def _apply_axis_flip_to_didq(didq: np.ndarray, axis_flip: tuple[int, int, int]) -> np.ndarray:
    out = np.zeros_like(didq)
    for a in range(3):
        for b in range(a + 1):
            idx = IDX_MAP[(a, b)]
            out[idx, :] = axis_flip[a] * axis_flip[b] * didq[idx, :]
    return out


def align_model_to_gaussian(log_path: str | Path, fchk_path: str | Path):
    didq_log_print = _parse_didq_block(log_path)
    force = parse_gaussian_anharmonic_force_data(log_path)
    harmonic = parse_gaussian_harmonic_data(log_path)
    fchk = parse_gaussian_fchk_harmonic_data(fchk_path)
    coords_ang = fchk.coords_bohr * 0.529177210903
    if harmonic.harmonic_to_anharmonic is None:
        raise ValueError("Gaussian harmonic-to-anharmonic equivalency table is required.")
    a_to_h = np.asarray(harmonic.harmonic_to_anharmonic, dtype=int)
    h_to_a = np.empty_like(a_to_h)
    for a_idx, h_idx in enumerate(a_to_h):
        h_to_a[h_idx] = a_idx
    didq_log = didq_log_print[:, a_to_h]

    best = None
    best_model = None
    for rep in ("I", "II", "III"):
        model = harmonic_inertia_model_from_geometry_hessian(
            fchk.masses_amu,
            coords_ang,
            fchk.cartesian_force_constants,
            representation=rep,
        )
        # The aligned working order is the Gaussian printed/anharmonic (A) order.
        # Reorder the harmonic model directly to the harmonic indices requested by
        # the Gaussian H<->A equivalency table, then compare to the printed dIdQ block.
        order = np.array(a_to_h, dtype=int)
        _apply_mode_reorder_and_signs(model, order, np.ones_like(order, dtype=float))
        didq0 = _didq_from_model(model)

        mode_signs = np.ones(didq0.shape[1], dtype=float)
        for k in range(didq0.shape[1]):
            if np.dot(didq0[:, k], didq_log[:, k]) < 0.0:
                mode_signs[k] = -1.0
        _apply_mode_reorder_and_signs(model, np.arange(didq0.shape[1], dtype=int), mode_signs)
        didq1 = _didq_from_model(model)

        for axis_flip in itertools.product((1, -1), repeat=3):
            didq2 = _apply_axis_flip_to_didq(didq1, axis_flip)
            err = float(np.linalg.norm(didq2 - didq_log))
            if best is None or err < best.didq_error:
                best = AlignmentResult(
                    representation=rep,
                    axis_flip=tuple(int(x) for x in axis_flip),
                    printed_to_harmonic=[int(x) + 1 for x in a_to_h],
                    harmonic_to_printed=[int(x) + 1 for x in h_to_a],
                    mode_signs=[int(x) for x in mode_signs],
                    didq_error=err,
                )
                best_model = model

    assert best is not None and best_model is not None
    return best_model, best, didq_log_print


def _axis_flipped_coords_and_disp(model, mode_index_1based: int, axis_flip: tuple[int, int, int]):
    idx = mode_index_1based - 1
    coords = np.asarray(model.coords_pa_ang, dtype=float).copy()
    masses = np.asarray(model.masses_amu, dtype=float)
    modes = np.asarray(model.vib_vecs_mw_pa, dtype=float).reshape(len(masses), 3, -1)
    disp = modes[:, :, idx] / np.sqrt(masses)[:, None]
    for a, s in enumerate(axis_flip):
        coords[:, a] *= s
        disp[:, a] *= s
    return masses, coords, disp


def best_axis_flip_for_mode(model, didq_log_print: np.ndarray, printed_mode_1based: int) -> tuple[tuple[int, int, int], float]:
    ref = didq_log_print[:, printed_mode_1based - 1]
    didq0 = _didq_from_model(model)[:, printed_mode_1based - 1]
    best_flip = (1, 1, 1)
    best_err = float("inf")
    for flip in itertools.product((1, -1), repeat=3):
        d = _apply_axis_flip_to_didq(didq0[:, None], flip)[:, 0]
        err = float(np.linalg.norm(d - ref))
        if err < best_err:
            best_flip = tuple(int(x) for x in flip)
            best_err = err
    return best_flip, best_err


def inertia_tensor(masses_amu: np.ndarray, coords_ang: np.ndarray) -> np.ndarray:
    x = coords_ang[:, 0]
    y = coords_ang[:, 1]
    z = coords_ang[:, 2]
    ixx = np.sum(masses_amu * (y * y + z * z))
    iyy = np.sum(masses_amu * (x * x + z * z))
    izz = np.sum(masses_amu * (x * x + y * y))
    ixy = -np.sum(masses_amu * x * y)
    ixz = -np.sum(masses_amu * x * z)
    iyz = -np.sum(masses_amu * y * z)
    return np.array([[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]], dtype=float)


def inertia_components_along_mode(model, alignment: AlignmentResult, mode_index_1based: int, q_values: np.ndarray) -> dict[str, np.ndarray]:
    masses, coords0, disp = _axis_flipped_coords_and_disp(model, mode_index_1based, alignment.axis_flip)
    comps = {label: np.zeros_like(q_values, dtype=float) for label in PAIR_LABELS}
    for iq, q in enumerate(q_values):
        tens = inertia_tensor(masses, coords0 + q * disp)
        comps["Ixx"][iq] = tens[0, 0]
        comps["Ixy"][iq] = tens[0, 1]
        comps["Iyy"][iq] = tens[1, 1]
        comps["Ixz"][iq] = tens[0, 2]
        comps["Iyz"][iq] = tens[1, 2]
        comps["Izz"][iq] = tens[2, 2]
    return comps


def rotational_constants_from_inertia_components(components: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for label in ("Ixx", "Iyy", "Izz"):
        inv = 1.0 / np.asarray(components[label], dtype=float)
        out[label.replace("I", "1/I")] = inv
        out[label.replace("I", "B_")] = B_MHZ_FROM_I_AMU_A2 * inv
    return out


def mode_potential(force, mode_index_1based: int, q_values: np.ndarray) -> dict[str, np.ndarray | float]:
    idx = mode_index_1based - 1
    omega = float(force.frequencies_cm[idx])
    phi3 = float(force.phi3_reduced_cm[idx, idx, idx])
    phi4 = float(force.phi4_reduced_cm[idx, idx, idx, idx])
    v = 0.5 * omega * omega * q_values * q_values
    v += (phi3 / 6.0) * q_values**3
    v += (phi4 / 24.0) * q_values**4
    return {"omega_cm": omega, "phi3_cm": phi3, "phi4_cm": phi4, "V_cm": v}


def finite_difference_didq(model, alignment: AlignmentResult, mode_index_1based: int, h: float = 1.0e-4) -> np.ndarray:
    masses, coords0, disp = _axis_flipped_coords_and_disp(model, mode_index_1based, alignment.axis_flip)
    plus = inertia_tensor(masses, coords0 + h * disp)
    minus = inertia_tensor(masses, coords0 - h * disp)
    d = (plus - minus) / (2.0 * h)
    return np.array([d[0, 0], d[0, 1], d[1, 1], d[0, 2], d[1, 2], d[2, 2]], dtype=float)


def _q_operator_matrix(omega_cm: float, basis_size: int) -> np.ndarray:
    q = np.zeros((basis_size, basis_size), dtype=float)
    scale = 1.0 / math.sqrt(2.0 * omega_cm)
    for n in range(basis_size - 1):
        val = math.sqrt(n + 1.0) * scale
        q[n, n + 1] = val
        q[n + 1, n] = val
    return q


def _operator_from_polynomial(qmat: np.ndarray, coeffs_low_to_high: np.ndarray) -> np.ndarray:
    basis_size = qmat.shape[0]
    out = np.zeros_like(qmat)
    qpow = np.eye(basis_size)
    for c in coeffs_low_to_high:
        out += c * qpow
        qpow = qpow @ qmat
    return out


def _fit_polynomial(q_values: np.ndarray, values: np.ndarray, degree: int) -> np.ndarray:
    return np.polynomial.polynomial.polyfit(q_values, values, deg=degree)


def variational_vs_perturbative_average(
    omega_cm: float,
    phi3_cm: float,
    phi4_cm: float,
    q_values: np.ndarray,
    observable_values: np.ndarray,
    *,
    basis_size: int = 32,
    poly_degree: int = 10,
) -> dict[str, float]:
    """Compare variational and second-order perturbative averages in a 1D HO basis."""
    qmat = _q_operator_matrix(omega_cm, basis_size)
    h0 = np.diag(omega_cm * (np.arange(basis_size, dtype=float) + 0.5))
    q2 = qmat @ qmat
    q3 = q2 @ qmat
    q4 = q2 @ q2
    vmat = (phi3_cm / 6.0) * q3 + (phi4_cm / 24.0) * q4
    h = h0 + vmat

    obs_coeff = _fit_polynomial(q_values, observable_values, poly_degree)
    omat = _operator_from_polynomial(qmat, obs_coeff)

    evals, evecs = np.linalg.eigh(h)
    psi0 = evecs[:, 0]
    var_avg = float(np.real_if_close(psi0.conj() @ (omat @ psi0)))

    e0 = h0[0, 0]
    c1 = np.zeros(basis_size, dtype=float)
    for n in range(1, basis_size):
        c1[n] = vmat[n, 0] / (e0 - h0[n, n])

    c2 = np.zeros(basis_size, dtype=float)
    for n in range(1, basis_size):
        acc = 0.0
        for m in range(1, basis_size):
            acc += vmat[n, m] * c1[m] / (e0 - h0[n, n])
        acc -= vmat[0, 0] * c1[n] / (e0 - h0[n, n])
        c2[n] = acc

    o00 = float(omat[0, 0])
    pert = o00
    pert += 2.0 * float(np.dot(omat[0, 1:], c1[1:]))
    pert += float(c1[1:] @ (omat[1:, 1:] @ c1[1:]))
    pert += 2.0 * float(np.dot(omat[0, 1:], c2[1:]))
    pert -= o00 * float(np.dot(c1[1:], c1[1:]))

    harm = o00
    return {
        "harmonic": harm,
        "perturbative": pert,
        "variational": var_avg,
        "delta_var_minus_pert": var_avg - pert,
        "delta_var_minus_harm": var_avg - harm,
    }


def _aligned_phi3(force) -> np.ndarray:
    """The aligned working order is the Gaussian printed/anharmonic (A) order."""
    return np.asarray(force.phi3_reduced_cm, dtype=float).copy()


def main() -> None:
    ap = argparse.ArgumentParser(description="Probe one Gaussian normal mode in the aligned dIdQ convention.")
    ap.add_argument("--log", required=True, help="Gaussian anharmonic log")
    ap.add_argument("--fchk", required=True, help="Gaussian fchk")
    ap.add_argument("--mode", required=True, type=int, help="1-based mode index in Gaussian printed/anharmonic order")
    ap.add_argument("--mode-order", choices=("printed", "harmonic"), default="printed", help="Interpret --mode in Gaussian printed (A) order or harmonic (H) order")
    ap.add_argument("--qmax", type=float, default=1.0, help="Maximum Q displacement")
    ap.add_argument("--nq", type=int, default=101, help="Number of Q points")
    ap.add_argument("--axis", choices=AXES, help="If provided, compare variational and perturbative B_tau(Q) averages")
    ap.add_argument("--basis-size", type=int, default=32, help="HO basis size for 1D variational/PT comparison")
    ap.add_argument("--poly-degree", type=int, default=10, help="Polynomial degree used to represent B_tau(Q)")
    ap.add_argument("--json", action="store_true", help="Emit JSON payload")
    args = ap.parse_args()

    model, alignment, didq_log = align_model_to_gaussian(args.log, args.fchk)
    if args.mode_order == "printed":
        printed_mode = args.mode
        harmonic_mode = alignment.printed_to_harmonic[args.mode - 1]
    else:
        harmonic_mode = args.mode
        printed_mode = alignment.harmonic_to_printed[args.mode - 1]
    internal_mode = printed_mode
    local_flip, local_flip_err = best_axis_flip_for_mode(model, didq_log, printed_mode)
    local_alignment = AlignmentResult(
        representation=alignment.representation,
        axis_flip=local_flip,
        printed_to_harmonic=alignment.printed_to_harmonic,
        harmonic_to_printed=alignment.harmonic_to_printed,
        mode_signs=alignment.mode_signs,
        didq_error=local_flip_err,
    )

    q = np.linspace(-args.qmax, args.qmax, args.nq)
    components = inertia_components_along_mode(model, local_alignment, internal_mode, q)
    invs = rotational_constants_from_inertia_components(components)
    force = parse_gaussian_anharmonic_force_data(args.log)
    pot = mode_potential(force, internal_mode, q)
    didq_fd = finite_difference_didq(model, local_alignment, internal_mode)
    didq_ref = didq_log[:, printed_mode - 1]

    payload = {
        "alignment": asdict(alignment),
        "local_alignment": asdict(local_alignment),
        "mode": args.mode,
        "mode_order": args.mode_order,
        "printed_mode": printed_mode,
        "harmonic_mode": harmonic_mode,
        "internal_mode": internal_mode,
        "frequency_cm": float(model.vib_freq_cm[internal_mode - 1]),
        "didq_fd": didq_fd.tolist(),
        "didq_gaussian": didq_ref.tolist(),
        "didq_diff_norm": float(np.linalg.norm(didq_fd - didq_ref)),
        "q_grid": q.tolist(),
        "inertia": {k: v.tolist() for k, v in components.items()},
        "inverse_inertia": {k: v.tolist() for k, v in invs.items()},
        "potential": {
            "omega_cm": float(pot["omega_cm"]),
            "phi3_cm": float(pot["phi3_cm"]),
            "phi4_cm": float(pot["phi4_cm"]),
            "V_cm": pot["V_cm"].tolist(),
        },
        "note": (
            "The Q coordinate and I(Q) convention are aligned to Gaussian dIdQ. "
            "The 1D potential uses Gaussian reduced force constants as printed. "
            "The 1D variational/PT comparison uses a harmonic-oscillator basis "
            "and a polynomial representation of B_tau(Q)."
        ),
    }

    if args.axis is not None:
        axis_map_inv = {"A": "1/Ixx", "B": "1/Iyy", "C": "1/Izz"}
        axis_map_B = {"A": "B_xx", "B": "B_yy", "C": "B_zz"}
        key_inv = axis_map_inv[args.axis]
        key_B = axis_map_B[args.axis]
        comp_inv = variational_vs_perturbative_average(
            float(pot["omega_cm"]),
            float(pot["phi3_cm"]),
            float(pot["phi4_cm"]),
            q,
            np.asarray(invs[key_inv], dtype=float),
            basis_size=args.basis_size,
            poly_degree=args.poly_degree,
        )
        comp_B = variational_vs_perturbative_average(
            float(pot["omega_cm"]),
            float(pot["phi3_cm"]),
            float(pot["phi4_cm"]),
            q,
            np.asarray(invs[key_B], dtype=float),
            basis_size=args.basis_size,
            poly_degree=args.poly_degree,
        )

        force_aligned = parse_gaussian_anharmonic_force_data(args.log)
        phi3_aligned = _aligned_phi3(force_aligned)
        alpha_parts = alpha_matrix_from_harmonic_and_cubic_cm(model, phi3_aligned)
        abc_axes = ("A", "B", "C")
        axis_idx = abc_axes.index(args.axis)
        mode_idx = internal_mode - 1
        alpha_total_cm = float(alpha_parts["alpha_total_cm_abc"][mode_idx, axis_idx])
        alpha_cor_cm = float(alpha_parts["alpha_coriolis_cm_abc"][mode_idx, axis_idx])
        alpha_noncor_cm = alpha_total_cm - alpha_cor_cm
        alpha_total_mhz = alpha_total_cm * 29979.2458
        alpha_cor_mhz = alpha_cor_cm * 29979.2458
        alpha_noncor_mhz = alpha_noncor_cm * 29979.2458
        delta_b_vpt2_noncor = -0.5 * alpha_noncor_mhz
        delta_b_var = comp_B["variational"] - comp_B["harmonic"]

        payload["axis_comparison"] = {
            "axis": args.axis,
            "observable_inverse_inertia": key_inv,
            "observable_B": key_B,
            "inverse_inertia": comp_inv,
            "B_mhz": comp_B,
            "gaussian_alpha_total_mhz": alpha_total_mhz,
            "gaussian_alpha_coriolis_mhz": alpha_cor_mhz,
            "gaussian_alpha_noncoriolis_mhz": alpha_noncor_mhz,
            "gaussian_delta_b0_noncoriolis_mhz": delta_b_vpt2_noncor,
            "variational_delta_b0_mhz": delta_b_var,
            "delta_var_minus_vpt2_noncor_mhz": delta_b_var - delta_b_vpt2_noncor,
        }

    if args.json:
        print(json.dumps(payload, indent=2))
        return

    print(
        f"Alignment: rep={alignment.representation}, global_flip={alignment.axis_flip}, global_err={alignment.didq_error:.6g}, "
        f"local_flip={local_alignment.axis_flip}, local_err={local_alignment.didq_error:.6g}"
    )
    print(
        f"Mode {args.mode} ({args.mode_order}; printed={printed_mode}, harmonic={harmonic_mode}): "
        f"freq={payload['frequency_cm']:.6f} cm^-1"
    )
    print("Gaussian dIdQ:", " ".join(f"{x: .6f}" for x in didq_ref))
    print("FD dIdQ      :", " ".join(f"{x: .6f}" for x in didq_fd))
    print(f"dIdQ diff norm: {payload['didq_diff_norm']:.6g}")
    print(
        "Potential:",
        f"omega={payload['potential']['omega_cm']:.6f} cm^-1,",
        f"phi3={payload['potential']['phi3_cm']:.6f} cm^-1,",
        f"phi4={payload['potential']['phi4_cm']:.6f} cm^-1",
    )
    if "axis_comparison" in payload:
        comp = payload["axis_comparison"]
        inv = comp["inverse_inertia"]
        bmhz = comp["B_mhz"]
        print(
            f"{comp['axis']} axis ({comp['observable_inverse_inertia']}):",
            f"harm={inv['harmonic']:.12f},",
            f"pert={inv['perturbative']:.12f},",
            f"var={inv['variational']:.12f},",
            f"var-pert={inv['delta_var_minus_pert']:.12f}",
        )
        print(
            f"{comp['axis']} axis ({comp['observable_B']}) [MHz]:",
            f"harm={bmhz['harmonic']:.6f},",
            f"pert={bmhz['perturbative']:.6f},",
            f"var={bmhz['variational']:.6f},",
            f"var-pert={bmhz['delta_var_minus_pert']:.6f}",
        )
        print(
            "Gaussian single-mode [MHz]:",
            f"alpha_total={comp['gaussian_alpha_total_mhz']:.6f},",
            f"alpha_cor={comp['gaussian_alpha_coriolis_mhz']:.6f},",
            f"alpha_noncor={comp['gaussian_alpha_noncoriolis_mhz']:.6f},",
            f"deltaB0_noncor={comp['gaussian_delta_b0_noncoriolis_mhz']:.6f}",
        )
        print(
            "1D comparison [MHz]:",
            f"deltaB0_var={comp['variational_delta_b0_mhz']:.6f},",
            f"var - VPT2(noncor)={comp['delta_var_minus_vpt2_noncor_mhz']:.6f}",
        )
    print("Grid written only to stdout/JSON; use --json for downstream analysis.")


if __name__ == "__main__":
    main()
