"""Point-group irrep labels used by frozen GIC definitions."""

from __future__ import annotations

import re

import numpy as np


_IRREP_ORDER_BY_POINT_GROUP = {
    "C1": ("A",),
    "CS": ("A'", "A''"),
    "CI": ("Ag", "Au"),
    "C2": ("A", "B"),
    "C2V": ("A1", "A2", "B1", "B2"),
    "D2": ("A", "B1", "B2", "B3"),
    "C2H": ("Ag", "Bg", "Au", "Bu"),
    "D2H": ("Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"),
    "C3V": ("A1", "A2", "E"),
    "C4V": ("A1", "A2", "B1", "B2", "E"),
    "C6V": ("A1", "A2", "B1", "B2", "E1", "E2"),
    "D3H": ("A1'", "A2'", "E'", "A1''", "A2''", "E''"),
    "D4H": (
        "A1g",
        "A2g",
        "B1g",
        "B2g",
        "Eg",
        "A1u",
        "A2u",
        "B1u",
        "B2u",
        "Eu",
    ),
    "TD": ("A1", "A2", "E", "T1", "T2"),
    "OH": ("A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu", "T1u", "T2u"),
    "IH": ("Ag", "T1g", "T2g", "Gg", "Hg", "Au", "T1u", "T2u", "Gu", "Hu"),
}


def normalized_point_group(point_group: str | None) -> str:
    text = (point_group or "C1").strip()
    return text or "C1"


def irrep_sequence(point_group: str | None) -> tuple[str, ...]:
    """Return a deterministic irrep order with the totally symmetric irrep first."""
    group = normalized_point_group(point_group)
    exact = _IRREP_ORDER_BY_POINT_GROUP.get(group.upper())
    if exact:
        return exact

    match = re.fullmatch(r"([CD])(\d+)([VHD]?)", group.upper())
    if not match:
        return ("A",)
    family, n_text, suffix = match.groups()
    n = int(n_text)
    if family == "C":
        if suffix == "V":
            base = ["A1", "A2"]
        elif suffix == "H":
            if n == 2:
                return _IRREP_ORDER_BY_POINT_GROUP["C2H"]
            if n % 2 == 1:
                base = ["A'", "A''"]
                for idx in range(1, (n + 1) // 2):
                    base.extend([f"{_e_name(idx, n)}'", f"{_e_name(idx, n)}''"])
                return tuple(dict.fromkeys(base))
            base = ["Ag", "Au"]
        else:
            base = ["A"]
        if n % 2 == 0:
            base.extend(["B1", "B2"] if suffix == "V" else ["B"])
        base.extend(f"E{idx}" for idx in range(1, max(1, n // 2)))
        return tuple(dict.fromkeys(base))

    base = ["A1", "A2"]
    if n % 2 == 0:
        base.extend(["B1", "B2"])
    base.extend(f"E{idx}" for idx in range(1, max(1, n // 2)))
    if suffix == "H" and n % 2 == 0:
        return tuple(
            dict.fromkeys(
                [f"{name}g" for name in base if not name.startswith("E")]
                + [name + "g" for name in base if name.startswith("E")]
                + [f"{name}u" for name in base if not name.startswith("E")]
                + [name + "u" for name in base if name.startswith("E")]
            )
        )
    if suffix == "H":
        return tuple(
            dict.fromkeys(
                [f"{name}'" for name in base if not name.startswith("E")]
                + [name + "'" for name in base if name.startswith("E")]
                + [f"{name}''" for name in base if not name.startswith("E")]
                + [name + "''" for name in base if name.startswith("E")]
            )
        )
    return tuple(dict.fromkeys(base))


def total_symmetric_irrep(point_group: str | None) -> str:
    return irrep_sequence(point_group)[0]


def non_total_irrep_sequence(point_group: str | None) -> tuple[str, ...]:
    irreps = irrep_sequence(point_group)
    return tuple(irrep for irrep in irreps if irrep != irreps[0])


def is_total_symmetric_irrep(point_group: str | None, irrep: str | None) -> bool:
    return (irrep or "").strip() == total_symmetric_irrep(point_group)


def irrep_name_prefix(irrep: str | None) -> str:
    """Convert an irrep label to a compact coordinate-name prefix."""
    text = (irrep or "X").strip() or "X"
    return text.replace("'", "p").replace('"', "pp")


def irrep_characters_for_operations(
    operation_labels: tuple[str, ...] | list[str],
    point_group: str | None,
) -> tuple[tuple[str, tuple[float, ...]], ...]:
    labels = tuple(_canonical_operation_label(label) for label in operation_labels)
    if labels == ("E",):
        return (("A", (1.0,)),)
    group = normalized_point_group(point_group)
    group_key = group.upper()
    if group_key == "CS":
        return (
            ("A'", tuple(1.0 for _label in labels)),
            (
                "A''",
                tuple(-1.0 if label.startswith("sigma") else 1.0 for label in labels),
            ),
        )
    if group_key == "CI":
        return (
            ("Ag", tuple(1.0 for _label in labels)),
            ("Au", tuple(-1.0 if label == "i" else 1.0 for label in labels)),
        )
    if group_key == "C2":
        return (
            ("A", tuple(1.0 for _label in labels)),
            ("B", tuple(-1.0 if label.startswith("C2") else 1.0 for label in labels)),
        )
    if group_key == "C2H":
        return _c2h_characters(labels)
    if group_key == "C2V":
        sigma_labels = [label for label in labels if label.startswith("sigma")]
        preferred = ("sigma_xz", "sigma_yz", "sigma_xy")
        first_sigma = next((label for label in preferred if label in sigma_labels), None)
        second_sigma = next(
            (label for label in sorted(sigma_labels) if label != first_sigma),
            None,
        )
        if first_sigma and second_sigma:
            table = {
                "E": (1.0, 1.0, 1.0, 1.0),
                "C2x": (1.0, 1.0, -1.0, -1.0),
                "C2y": (1.0, 1.0, -1.0, -1.0),
                "C2z": (1.0, 1.0, -1.0, -1.0),
                "C2_perp": (1.0, 1.0, -1.0, -1.0),
                first_sigma: (1.0, -1.0, 1.0, -1.0),
                second_sigma: (1.0, -1.0, -1.0, 1.0),
            }
            rows = np.array(
                [
                    table.get(
                        "C2_perp" if label.startswith("C2_xy") else label,
                        (1.0, 1.0, 1.0, 1.0),
                    )
                    for label in labels
                ]
            )
            return tuple(
                (name, tuple(float(value) for value in rows[:, idx]))
                for idx, name in enumerate(("A1", "A2", "B1", "B2"))
            )
    if group_key == "D2":
        return _d2_characters(labels)
    if group_key == "D2H":
        return _d2h_characters(labels)
    generic = _generic_family_characters(labels, group_key)
    if generic:
        return generic
    return ()


def _generic_family_characters(
    labels: tuple[str, ...],
    group_key: str,
) -> tuple[tuple[str, tuple[float, ...]], ...]:
    match = re.fullmatch(r"([CD])(\d+)([VHD]?)", group_key)
    if not match:
        return ()
    family, n_text, suffix = match.groups()
    n = int(n_text)
    if n < 2:
        return ()
    if family == "C" and suffix == "":
        return _cn_characters(labels, n)
    if family == "C" and suffix == "V":
        return _cnv_characters(labels, n)
    if family == "C" and suffix == "H":
        return _cnh_characters(labels, n) if n % 2 == 1 else ()
    if family == "D" and suffix == "":
        return _dn_characters(labels, n)
    if family == "D" and suffix == "H":
        return _dnh_characters(labels, n)
    return ()


def _cn_characters(
    labels: tuple[str, ...],
    n: int,
) -> tuple[tuple[str, tuple[float, ...]], ...]:
    rows: list[tuple[str, tuple[float, ...]]] = []
    rows.append(("A", tuple(_rotation_character(label, n, 0) for label in labels)))
    if n % 2 == 0:
        rows.append(("B", tuple(_rotation_character(label, n, n // 2) for label in labels)))
    for order in range(1, (n + 1) // 2):
        if n % 2 == 0 and order == n // 2:
            continue
        rows.append(
            (
                _e_name(order, n),
                tuple(2.0 * _rotation_character(label, n, order) for label in labels),
            )
        )
    return tuple(rows)


def _cnv_characters(
    labels: tuple[str, ...],
    n: int,
) -> tuple[tuple[str, tuple[float, ...]], ...]:
    rows: list[tuple[str, tuple[float, ...]]] = [
        ("A1", tuple(_cnv_one_dim_character(label, n, 0, 1) for label in labels)),
        ("A2", tuple(_cnv_one_dim_character(label, n, 0, -1) for label in labels)),
    ]
    if n % 2 == 0:
        rows.extend(
            (
                (
                    "B1",
                    tuple(_cnv_one_dim_character(label, n, n // 2, 1) for label in labels),
                ),
                (
                    "B2",
                    tuple(_cnv_one_dim_character(label, n, n // 2, -1) for label in labels),
                ),
            )
        )
    for order in range(1, (n + 1) // 2):
        if n % 2 == 0 and order == n // 2:
            continue
        rows.append(
            (
                _e_name(order, n),
                tuple(_cnv_e_character(label, n, order) for label in labels),
            )
        )
    return tuple(rows)


def _cnh_characters(
    labels: tuple[str, ...],
    n: int,
) -> tuple[tuple[str, tuple[float, ...]], ...]:
    rows: list[tuple[str, tuple[float, ...]]] = []
    base = [("A", 0)]
    if n % 2 == 0:
        base.append(("B", n // 2))
    for name, order in base:
        rows.append((name + "'", tuple(_cnh_character(label, n, order, 1) for label in labels)))
        rows.append((name + "''", tuple(_cnh_character(label, n, order, -1) for label in labels)))
    for order in range(1, (n + 1) // 2):
        if n % 2 == 0 and order == n // 2:
            continue
        rows.append(
            (
                _e_name(order, n) + "'",
                tuple(2.0 * _cnh_character(label, n, order, 1) for label in labels),
            )
        )
        rows.append(
            (
                _e_name(order, n) + "''",
                tuple(2.0 * _cnh_character(label, n, order, -1) for label in labels),
            )
        )
    return tuple(rows)


def _dn_characters(
    labels: tuple[str, ...],
    n: int,
) -> tuple[tuple[str, tuple[float, ...]], ...]:
    if n == 2:
        return _d2_characters(labels)
    rows: list[tuple[str, tuple[float, ...]]] = [
        ("A1", tuple(_dn_one_dim_character(label, n, 0, 1) for label in labels)),
        ("A2", tuple(_dn_one_dim_character(label, n, 0, -1) for label in labels)),
    ]
    if n % 2 == 0:
        rows.extend(
            (
                (
                    "B1",
                    tuple(_dn_one_dim_character(label, n, n // 2, 1) for label in labels),
                ),
                (
                    "B2",
                    tuple(_dn_one_dim_character(label, n, n // 2, -1) for label in labels),
                ),
            )
        )
    for order in range(1, (n + 1) // 2):
        if n % 2 == 0 and order == n // 2:
            continue
        rows.append(
            (
                _e_name(order, n),
                tuple(_dn_e_character(label, n, order) for label in labels),
            )
        )
    return tuple(rows)


def _rotation_character(label: str, n: int, order: int) -> float:
    power = _principal_rotation_power(label, n)
    if power is None:
        return 0.0
    return float(np.cos(2.0 * np.pi * order * power / n))


def _cnv_one_dim_character(label: str, n: int, order: int, reflection_sign: int) -> float:
    if _is_vertical_reflection(label):
        mirror_index = _vertical_reflection_index(label)
        phase = 1.0 if order == 0 else (-1.0 if (mirror_index % 2) else 1.0)
        return float(reflection_sign) * phase
    return _rotation_character(label, n, order)


def _cnv_e_character(label: str, n: int, order: int) -> float:
    if _is_vertical_reflection(label):
        return 0.0
    return 2.0 * _rotation_character(label, n, order)


def _cnh_character(label: str, n: int, order: int, reflection_sign: int) -> float:
    reflected = _is_horizontal_reflected(label, n)
    power = _cnh_rotation_power(label, n)
    if power is None:
        return 0.0
    value = float(np.cos(2.0 * np.pi * order * power / n))
    return value * (float(reflection_sign) if reflected else 1.0)


def _dn_one_dim_character(label: str, n: int, order: int, c2_sign: int) -> float:
    if _is_c2_prime(label):
        mirror_index = _c2_prime_index(label)
        phase = 1.0 if order == 0 else (-1.0 if (mirror_index % 2) else 1.0)
        return float(c2_sign) * phase
    return _rotation_character(label, n, order)


def _dn_e_character(label: str, n: int, order: int) -> float:
    if _is_c2_prime(label):
        return 0.0
    return 2.0 * _rotation_character(label, n, order)


def _dnh_characters(
    labels: tuple[str, ...],
    n: int,
) -> tuple[tuple[str, tuple[float, ...]], ...]:
    if n == 2:
        return _d2h_characters(labels)
    base_names = ["A1", "A2"]
    if n % 2 == 0:
        base_names.extend(["B1", "B2"])
    for order in range(1, (n + 1) // 2):
        if n % 2 == 0 and order == n // 2:
            continue
        base_names.append(_e_name(order, n))

    if n % 2 == 0:
        return tuple(
            (name + suffix, tuple(_dnh_value(label, n, name, parity) for label in labels))
            for suffix, parity in (("g", 1), ("u", -1))
            for name in base_names
        )
    return tuple(
        (name + suffix, tuple(_dnh_value(label, n, name, parity) for label in labels))
        for suffix, parity in (("'", 1), ("''", -1))
        for name in base_names
    )


def _dnh_value(label: str, n: int, irrep: str, reflection_sign: int) -> float:
    underlying, reflected = _dnh_underlying_operation(label, n)
    value = _dn_irrep_value(underlying, n, irrep)
    return value * (float(reflection_sign) if reflected else 1.0)


def _dnh_underlying_operation(label: str, n: int) -> tuple[str, bool]:
    if label == "sigma_xy":
        return "E", True
    if label == "i" and n % 2 == 0:
        return (f"C{n}z^{n // 2}", True)
    if _is_vertical_reflection(label):
        return label.replace("sigma_v", "C2_xy"), True
    match = re.fullmatch(r"sigma_h\*C(\d+)z\^(\d+)", label)
    if match:
        order, power = match.groups()
        return f"C{order}z^{power}", True
    return label, False


def _dn_irrep_value(label: str, n: int, irrep: str) -> float:
    if irrep == "A1":
        return _dn_one_dim_character(label, n, 0, 1)
    if irrep == "A2":
        return _dn_one_dim_character(label, n, 0, -1)
    if irrep == "B1" and n % 2 == 0:
        return _dn_one_dim_character(label, n, n // 2, 1)
    if irrep == "B2" and n % 2 == 0:
        return _dn_one_dim_character(label, n, n // 2, -1)
    match = re.fullmatch(r"E(\d*)", irrep)
    if match:
        order = int(match.group(1) or "1")
        return _dn_e_character(label, n, order)
    return 0.0


def _d2_characters(labels: tuple[str, ...]) -> tuple[tuple[str, tuple[float, ...]], ...]:
    values = {
        "A": {"C2x": 1.0, "C2y": 1.0, "C2z": 1.0},
        "B1": {"C2x": -1.0, "C2y": -1.0, "C2z": 1.0},
        "B2": {"C2x": -1.0, "C2y": 1.0, "C2z": -1.0},
        "B3": {"C2x": 1.0, "C2y": -1.0, "C2z": -1.0},
    }
    return tuple(
        (name, tuple(chars.get(label, 1.0) for label in labels))
        for name, chars in values.items()
    )


def _c2h_characters(labels: tuple[str, ...]) -> tuple[tuple[str, tuple[float, ...]], ...]:
    values = {
        "Ag": {"E": 1.0, "C2z": 1.0, "i": 1.0, "sigma_xy": 1.0},
        "Bg": {"E": 1.0, "C2z": -1.0, "i": 1.0, "sigma_xy": -1.0},
        "Au": {"E": 1.0, "C2z": 1.0, "i": -1.0, "sigma_xy": -1.0},
        "Bu": {"E": 1.0, "C2z": -1.0, "i": -1.0, "sigma_xy": 1.0},
    }
    return tuple(
        (name, tuple(chars.get(label, 0.0) for label in labels))
        for name, chars in values.items()
    )


def _d2h_characters(labels: tuple[str, ...]) -> tuple[tuple[str, tuple[float, ...]], ...]:
    rows: list[tuple[str, tuple[float, ...]]] = []
    for name in ("A", "B1", "B2", "B3"):
        rows.append((name + "g", tuple(_d2h_value(name, label, 1) for label in labels)))
    for name in ("A", "B1", "B2", "B3"):
        rows.append((name + "u", tuple(_d2h_value(name, label, -1) for label in labels)))
    return tuple(rows)


def _d2h_value(name: str, label: str, inversion_sign: int) -> float:
    d2_values = {
        "A": {"E": 1.0, "C2z": 1.0, "C2y": 1.0, "C2x": 1.0},
        "B1": {"E": 1.0, "C2z": 1.0, "C2y": -1.0, "C2x": -1.0},
        "B2": {"E": 1.0, "C2z": -1.0, "C2y": 1.0, "C2x": -1.0},
        "B3": {"E": 1.0, "C2z": -1.0, "C2y": -1.0, "C2x": 1.0},
    }[name]
    if label in d2_values:
        return d2_values[label]
    if label == "i":
        return float(inversion_sign)
    sigma_to_rotation = {
        "sigma_xy": "C2z",
        "sigma_xz": "C2y",
        "sigma_yz": "C2x",
    }
    if label in sigma_to_rotation:
        return float(inversion_sign) * d2_values[sigma_to_rotation[label]]
    return 0.0


def _e_name(order: int, n: int) -> str:
    if order == 1 and n in {3, 4}:
        return "E"
    return f"E{order}"


def _principal_rotation_power(label: str, n: int) -> int | None:
    if label == "E":
        return 0
    match = re.fullmatch(r"C(\d+)([xyz])\^(\d+)", label)
    if match:
        order, axis, power = match.groups()
        if axis != "z":
            return None
        order_int = int(order)
        power_int = int(power)
        relative = n * power_int / order_int
        nearest = int(round(relative))
        if abs(relative - nearest) > 1.0e-8:
            return None
        return nearest % n
    match = re.fullmatch(r"C2([xyz])", label)
    if match:
        axis = match.group(1)
        if axis == "z" and n % 2 == 0:
            return n // 2
    return None


def _cnh_rotation_power(label: str, n: int) -> int | None:
    if label == "sigma_xy":
        return 0
    if label == "i" and n % 2 == 0:
        return n // 2
    match = re.fullmatch(r"sigma_h\*C(\d+)z\^(\d+)", label)
    if match:
        order, power = match.groups()
        relative = n * int(power) / int(order)
        nearest = int(round(relative))
        if abs(relative - nearest) > 1.0e-8:
            return None
        return nearest % n
    return _principal_rotation_power(label, n)


def _is_vertical_reflection(label: str) -> bool:
    return label in {"sigma_xz", "sigma_yz"} or label.startswith("sigma_v")


def _vertical_reflection_index(label: str) -> int:
    match = re.fullmatch(r"sigma_v_(\d+)_(\d+)", label)
    if match:
        return int(match.group(2))
    if label == "sigma_xz":
        return 0
    if label == "sigma_yz":
        return 1
    return 0


def _is_horizontal_reflected(label: str, n: int) -> bool:
    if label == "sigma_xy" or label.startswith("sigma_h"):
        return True
    return bool(label == "i" and n % 2 == 0)


def _is_c2_prime(label: str) -> bool:
    return label.startswith("C2_xy") or label in {"C2x", "C2y"}


def _c2_prime_index(label: str) -> int:
    match = re.fullmatch(r"C2_xy_(\d+)_(\d+)", label)
    if match:
        return int(match.group(2))
    if label == "C2x":
        return 0
    if label == "C2y":
        return 1
    return 0


def _canonical_operation_label(label: str) -> str:
    text = str(label)
    if text in {"E", "i"} or text.startswith("sigma"):
        return text
    match = re.match(r"C(\d+)([xyz])\^(\d+)", text)
    if match:
        order, axis, _power = match.groups()
        return f"C2{axis}" if int(order) == 2 else text
    return text
