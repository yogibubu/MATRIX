from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote, urljoin
from urllib.request import Request, urlopen
import csv
import math
import re

import numpy as np

from .contracts import VibrationalSection, read_vibrational_section


NIST_WEBBOOK_BASE = "https://webbook.nist.gov"
NIST_USER_AGENT = "ORACLE-MATRIX/0.1"
VIBRATIONAL_OBSERVABLES = ("IR", "RAMAN", "VCD", "ROA")


@dataclass(frozen=True)
class VibrationalPeak:
    mode: int
    frequency_cm1: float
    intensity: float
    observable: str
    source: str


@dataclass(frozen=True)
class VibrationalSpectrumOptions:
    fwhm_cm1: float = 10.0
    step_cm1: float = 1.0
    lineshape: str = "gaussian"
    normalize: bool = True
    x_min_cm1: float | None = None
    x_max_cm1: float | None = None


@dataclass(frozen=True)
class VibrationalSpectrum:
    observable: str
    source: str
    x_cm1: np.ndarray
    y: np.ndarray
    peaks: tuple[VibrationalPeak, ...]
    y_label: str


@dataclass(frozen=True)
class VibrationalSpectrumComparison:
    first: VibrationalSpectrum
    second: VibrationalSpectrum
    x_cm1: np.ndarray
    first_y: np.ndarray
    second_y: np.ndarray
    plotted_second_y: np.ndarray
    mirror_second: bool


@dataclass(frozen=True)
class NISTIRPoint:
    wavenumber_cm1: float
    value: float


@dataclass(frozen=True)
class NISTIRDownloadResult:
    status: str
    message: str
    identifier: str
    page_url: str
    jcamp_url: str = ""
    state: str = ""
    csv_path: Path | None = None
    points: tuple[NISTIRPoint, ...] = ()

    @property
    def needs_user_instruction(self) -> bool:
        return self.status in {"not_gas_phase", "not_found", "parse_error"}


def vibrational_peaks_from_section(
    section: VibrationalSection,
    *,
    observable: str = "IR",
    source: str = "harmonic",
) -> tuple[VibrationalPeak, ...]:
    obs = _normalize_observable(observable)
    src = _normalize_source(source)
    frequencies = _frequencies_for_source(section, src)
    intensities = _intensities_for_observable(section, obs, src)
    if not frequencies:
        raise ValueError(f"no {src} vibrational frequencies are available")
    if not intensities:
        raise ValueError(f"no {src} {obs} intensities are available")
    if len(intensities) != len(frequencies):
        raise ValueError(
            f"{src} {obs} intensity count ({len(intensities)}) does not match "
            f"frequency count ({len(frequencies)})"
        )
    return tuple(
        VibrationalPeak(
            mode=idx,
            frequency_cm1=float(freq),
            intensity=float(intensity),
            observable=obs,
            source=src,
        )
        for idx, (freq, intensity) in enumerate(zip(frequencies, intensities), start=1)
    )


def build_vibrational_spectrum(
    section: VibrationalSection,
    *,
    observable: str = "IR",
    source: str = "harmonic",
    options: VibrationalSpectrumOptions | None = None,
) -> VibrationalSpectrum:
    opts = options or VibrationalSpectrumOptions()
    peaks = vibrational_peaks_from_section(section, observable=observable, source=source)
    if opts.fwhm_cm1 <= 0.0:
        raise ValueError("fwhm_cm1 must be positive")
    if opts.step_cm1 <= 0.0:
        raise ValueError("step_cm1 must be positive")
    centers = np.asarray([peak.frequency_cm1 for peak in peaks], dtype=float)
    x_min = (
        float(opts.x_min_cm1)
        if opts.x_min_cm1 is not None
        else max(0.0, centers.min() - 5.0 * opts.fwhm_cm1)
    )
    x_max = (
        float(opts.x_max_cm1) if opts.x_max_cm1 is not None else centers.max() + 5.0 * opts.fwhm_cm1
    )
    if x_max <= x_min:
        raise ValueError("x_max_cm1 must be greater than x_min_cm1")
    x = np.arange(x_min, x_max + 0.5 * opts.step_cm1, opts.step_cm1, dtype=float)
    y = np.zeros_like(x)
    for peak in peaks:
        y += peak.intensity * _profile(x, peak.frequency_cm1, opts.fwhm_cm1, opts.lineshape)
    if opts.normalize:
        scale = float(np.max(np.abs(y))) if y.size else 0.0
        if scale > 0.0:
            y = y / scale
    return VibrationalSpectrum(
        observable=_normalize_observable(observable),
        source=_normalize_source(source),
        x_cm1=x,
        y=y,
        peaks=peaks,
        y_label="normalized intensity" if opts.normalize else _y_label(observable),
    )


def build_vibrational_spectrum_from_xyzin(
    xyzin: Path | str,
    *,
    observable: str = "IR",
    source: str = "harmonic",
    options: VibrationalSpectrumOptions | None = None,
) -> VibrationalSpectrum:
    return build_vibrational_spectrum(
        read_vibrational_section(Path(xyzin)),
        observable=observable,
        source=source,
        options=options,
    )


def write_vibrational_spectrum_csv(path: Path | str, spectrum: VibrationalSpectrum) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavenumber_cm-1", spectrum.y_label, "observable", "source"])
        for x, y in zip(spectrum.x_cm1, spectrum.y):
            writer.writerow(
                [f"{float(x):.8f}", f"{float(y):.12g}", spectrum.observable, spectrum.source]
            )
    return target


def write_vibrational_peak_csv(path: Path | str, peaks: tuple[VibrationalPeak, ...]) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mode", "frequency_cm-1", "intensity", "observable", "source"])
        for peak in peaks:
            writer.writerow(
                [
                    peak.mode,
                    f"{peak.frequency_cm1:.8f}",
                    f"{peak.intensity:.12g}",
                    peak.observable,
                    peak.source,
                ]
            )
    return target


def write_vibrational_spectrum_plot(path: Path | str, spectrum: VibrationalSpectrum) -> Path:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    ax.plot(spectrum.x_cm1, spectrum.y, color="#225ea8", linewidth=1.4)
    ax.axhline(0.0, color="#4d4d4d", linewidth=0.6)
    ax.set_xlabel("Wavenumber / cm$^{-1}$")
    ax.set_ylabel(spectrum.y_label)
    ax.set_title(f"{spectrum.source.capitalize()} {spectrum.observable} spectrum")
    ax.invert_xaxis()
    fig.savefig(target)
    plt.close(fig)
    return target


def write_vibrational_spectrum_outputs(
    xyzin: Path | str,
    *,
    csv_path: Path | str,
    plot_path: Path | str | None = None,
    peaks_path: Path | str | None = None,
    observable: str = "IR",
    source: str = "harmonic",
    options: VibrationalSpectrumOptions | None = None,
) -> VibrationalSpectrum:
    spectrum = build_vibrational_spectrum_from_xyzin(
        xyzin,
        observable=observable,
        source=source,
        options=options,
    )
    write_vibrational_spectrum_csv(csv_path, spectrum)
    if peaks_path is not None:
        write_vibrational_peak_csv(peaks_path, spectrum.peaks)
    if plot_path is not None:
        write_vibrational_spectrum_plot(plot_path, spectrum)
    return spectrum


def compare_vibrational_spectra(
    first: VibrationalSpectrum,
    second: VibrationalSpectrum,
    *,
    mirror_second: bool | None = None,
) -> VibrationalSpectrumComparison:
    if first.x_cm1.size == 0 or second.x_cm1.size == 0:
        raise ValueError("cannot compare empty spectra")
    auto_mirror = _should_mirror_second(first.observable, second.observable)
    mirror = auto_mirror if mirror_second is None else bool(mirror_second)
    x_min = max(float(np.min(first.x_cm1)), float(np.min(second.x_cm1)))
    x_max = min(float(np.max(first.x_cm1)), float(np.max(second.x_cm1)))
    if x_max <= x_min:
        raise ValueError("spectra do not overlap in wavenumber")
    step = min(_grid_step(first.x_cm1), _grid_step(second.x_cm1))
    x = np.arange(x_min, x_max + 0.5 * step, step, dtype=float)
    first_y = np.interp(x, first.x_cm1, first.y)
    second_y = np.interp(x, second.x_cm1, second.y)
    plotted_second = -second_y if mirror else second_y
    return VibrationalSpectrumComparison(
        first=first,
        second=second,
        x_cm1=x,
        first_y=first_y,
        second_y=second_y,
        plotted_second_y=plotted_second,
        mirror_second=mirror,
    )


def write_vibrational_spectrum_comparison_csv(
    path: Path | str,
    comparison: VibrationalSpectrumComparison,
) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "wavenumber_cm-1",
                "first_y",
                "second_y",
                "plotted_second_y",
                "first_observable",
                "first_source",
                "second_observable",
                "second_source",
                "mirror_second",
            ]
        )
        for x, first_y, second_y, plotted_second_y in zip(
            comparison.x_cm1,
            comparison.first_y,
            comparison.second_y,
            comparison.plotted_second_y,
        ):
            writer.writerow(
                [
                    f"{float(x):.8f}",
                    f"{float(first_y):.12g}",
                    f"{float(second_y):.12g}",
                    f"{float(plotted_second_y):.12g}",
                    comparison.first.observable,
                    comparison.first.source,
                    comparison.second.observable,
                    comparison.second.source,
                    int(comparison.mirror_second),
                ]
            )
    return target


def write_vibrational_spectrum_comparison_plot(
    path: Path | str,
    comparison: VibrationalSpectrumComparison,
    *,
    first_label: str | None = None,
    second_label: str | None = None,
) -> Path:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    first_text = first_label or f"{comparison.first.source} {comparison.first.observable}"
    second_text = second_label or f"{comparison.second.source} {comparison.second.observable}"
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    ax.plot(comparison.x_cm1, comparison.first_y, color="#225ea8", linewidth=1.4, label=first_text)
    ax.plot(
        comparison.x_cm1,
        comparison.plotted_second_y,
        color="#c44e52",
        linewidth=1.4,
        label=second_text,
    )
    ax.axhline(0.0, color="#4d4d4d", linewidth=0.7)
    ax.set_xlabel("Wavenumber / cm$^{-1}$")
    ax.set_ylabel("normalized intensity")
    ax.set_title("Vibrational spectrum comparison")
    ax.legend(frameon=False)
    ax.invert_xaxis()
    fig.savefig(target)
    plt.close(fig)
    return target


def write_vibrational_spectrum_comparison_outputs(
    xyzin: Path | str,
    *,
    csv_path: Path | str,
    plot_path: Path | str | None = None,
    observable: str = "IR",
    first_source: str = "harmonic",
    second_source: str = "anharmonic",
    options: VibrationalSpectrumOptions | None = None,
    mirror_second: bool | None = None,
) -> VibrationalSpectrumComparison:
    first = build_vibrational_spectrum_from_xyzin(
        xyzin,
        observable=observable,
        source=first_source,
        options=options,
    )
    second = build_vibrational_spectrum_from_xyzin(
        xyzin,
        observable=observable,
        source=second_source,
        options=options,
    )
    comparison = compare_vibrational_spectra(first, second, mirror_second=mirror_second)
    write_vibrational_spectrum_comparison_csv(csv_path, comparison)
    if plot_path is not None:
        write_vibrational_spectrum_comparison_plot(plot_path, comparison)
    return comparison


def nist_ir_points_to_spectrum(
    points: tuple[NISTIRPoint, ...],
    *,
    source: str = "nist-gas-experiment",
    normalize: bool = True,
) -> VibrationalSpectrum:
    if not points:
        raise ValueError("NIST IR point list is empty")
    x = np.asarray([point.wavenumber_cm1 for point in points], dtype=float)
    transmittance = np.asarray([point.value for point in points], dtype=float)
    y = 1.0 - transmittance
    if normalize:
        scale = float(np.max(np.abs(y))) if y.size else 0.0
        if scale > 0.0:
            y = y / scale
    return VibrationalSpectrum(
        observable="IR",
        source=source,
        x_cm1=x,
        y=y,
        peaks=(),
        y_label="normalized absorbance" if normalize else "absorbance proxy",
    )


def fetch_nist_ir_gas_phase_csv(
    identifier: str,
    csv_path: Path | str,
    *,
    index: int = 1,
    timeout: float = 20.0,
) -> NISTIRDownloadResult:
    page_url = _nist_ir_page_url(identifier, index=index)
    page = _fetch_text(page_url, timeout=timeout)
    jcamp_url = _find_jcamp_url(page, page_url)
    if not jcamp_url:
        return NISTIRDownloadResult(
            status="not_found",
            message="NIST IR spectrum not found; user instruction is required",
            identifier=identifier,
            page_url=page_url,
        )
    jcamp = _fetch_text(jcamp_url, timeout=timeout, encoding="latin1")
    metadata = _jcamp_metadata(jcamp)
    state = metadata.get("STATE", "").strip()
    if not state.upper().startswith("GAS"):
        return NISTIRDownloadResult(
            status="not_gas_phase",
            message=f"NIST IR spectrum is not gas phase (state={state or 'unknown'})",
            identifier=identifier,
            page_url=page_url,
            jcamp_url=jcamp_url,
            state=state,
        )
    points = parse_nist_jcamp_ir_points(jcamp)
    if not points:
        return NISTIRDownloadResult(
            status="parse_error",
            message="NIST JCAMP spectrum could not be converted to numeric CSV",
            identifier=identifier,
            page_url=page_url,
            jcamp_url=jcamp_url,
            state=state,
        )
    target = write_nist_ir_csv(csv_path, points, state=state, source_url=jcamp_url)
    return NISTIRDownloadResult(
        status="downloaded",
        message=f"Downloaded NIST gas-phase IR spectrum ({len(points)} points)",
        identifier=identifier,
        page_url=page_url,
        jcamp_url=jcamp_url,
        state=state,
        csv_path=target,
        points=points,
    )


def parse_nist_jcamp_ir_points(text: str) -> tuple[NISTIRPoint, ...]:
    metadata = _jcamp_metadata(text)
    xfactor = float(metadata.get("XFACTOR", "1") or "1")
    yfactor = float(metadata.get("YFACTOR", "1") or "1")
    deltax = float(metadata.get("DELTAX", "0") or "0") * xfactor
    points: list[NISTIRPoint] = []
    in_xy = False
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("##"):
            key = line[2:].split("=", 1)[0].strip().upper()
            in_xy = key == "XYDATA"
            if key == "END":
                break
            continue
        if not in_xy:
            continue
        numbers = [float(item) for item in _NUMBER_RE.findall(line)]
        if len(numbers) < 2:
            continue
        x0 = numbers[0] * xfactor
        for offset, y_raw in enumerate(numbers[1:]):
            x = x0 + offset * deltax if deltax else x0
            points.append(NISTIRPoint(wavenumber_cm1=x, value=y_raw * yfactor))
    return tuple(points)


def write_nist_ir_csv(
    path: Path | str,
    points: tuple[NISTIRPoint, ...],
    *,
    state: str,
    source_url: str,
) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["# source", source_url])
        writer.writerow(["# state", state])
        writer.writerow(["wavenumber_cm-1", "transmittance"])
        for point in points:
            writer.writerow([f"{point.wavenumber_cm1:.8f}", f"{point.value:.12g}"])
    return target


def _profile(x: np.ndarray, center: float, fwhm: float, lineshape: str) -> np.ndarray:
    kind = lineshape.strip().lower()
    if kind == "gaussian":
        sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        return np.exp(-0.5 * ((x - center) / sigma) ** 2)
    if kind == "lorentzian":
        gamma = 0.5 * fwhm
        return (gamma * gamma) / ((x - center) ** 2 + gamma * gamma)
    raise ValueError("lineshape must be 'gaussian' or 'lorentzian'")


def _should_mirror_second(first_observable: str, second_observable: str) -> bool:
    signed = {"VCD", "ROA"}
    return (
        _normalize_observable(first_observable) not in signed
        and _normalize_observable(second_observable) not in signed
    )


def _grid_step(x: np.ndarray) -> float:
    if x.size < 2:
        return 1.0
    diffs = np.diff(np.sort(x))
    positive = diffs[diffs > 0.0]
    if positive.size == 0:
        return 1.0
    return float(np.min(positive))


def _frequencies_for_source(section: VibrationalSection, source: str) -> tuple[float, ...]:
    if source == "anharmonic":
        return section.anharmonic_frequencies_cm1
    return section.frequencies_cm1


def _intensities_for_observable(
    section: VibrationalSection,
    observable: str,
    source: str,
) -> tuple[float, ...]:
    if source == "anharmonic":
        anharmonic = {
            "IR": section.anharmonic_ir_intensities_km_mol,
            "RAMAN": section.anharmonic_raman_activities_A4_amu,
            "VCD": section.anharmonic_vcd_rot_strengths,
            "ROA": section.anharmonic_roa_intensities,
        }[observable]
        if anharmonic:
            return anharmonic
    return {
        "IR": section.ir_intensities_km_mol,
        "RAMAN": section.raman_activities_A4_amu,
        "VCD": section.vcd_rot_strengths,
        "ROA": section.roa_intensities,
    }[observable]


def _normalize_observable(value: str) -> str:
    observable = value.strip().upper()
    if observable not in VIBRATIONAL_OBSERVABLES:
        raise ValueError(f"unsupported vibrational observable: {value}")
    return observable


def _normalize_source(value: str) -> str:
    source = value.strip().lower()
    if source not in {"harmonic", "anharmonic"}:
        raise ValueError("source must be 'harmonic' or 'anharmonic'")
    return source


def _y_label(observable: str) -> str:
    return {
        "IR": "IR intensity / km mol$^{-1}$",
        "RAMAN": "Raman activity / A$^4$ amu$^{-1}$",
        "VCD": "VCD rotational strength",
        "ROA": "ROA intensity",
    }[_normalize_observable(observable)]


def _nist_ir_page_url(identifier: str, *, index: int) -> str:
    ident = identifier.strip()
    if re.fullmatch(r"C?\d+", ident, flags=re.IGNORECASE):
        nist_id = ident.upper() if ident.upper().startswith("C") else f"C{ident}"
        query = f"ID={quote(nist_id)}"
    elif re.fullmatch(r"\d{2,7}-\d{2}-\d", ident):
        query = f"ID={quote(ident)}"
    else:
        query = f"Name={quote(ident)}"
    return f"{NIST_WEBBOOK_BASE}/cgi/cbook.cgi?{query}&Index={int(index)}&Type=IR-SPEC"


def _fetch_text(url: str, *, timeout: float, encoding: str = "utf-8") -> str:
    request = Request(url, headers={"User-Agent": NIST_USER_AGENT})
    with urlopen(request, timeout=timeout) as response:
        return response.read().decode(encoding, errors="replace")


def _find_jcamp_url(page: str, page_url: str) -> str:
    match = re.search(r'href="([^"]*JCAMP[^"]*)"', page)
    if not match:
        return ""
    return urljoin(page_url, match.group(1).replace("&amp;", "&"))


def _jcamp_metadata(text: str) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line.startswith("##") or "=" not in line:
            continue
        key, value = line[2:].split("=", 1)
        metadata[key.strip().upper()] = value.strip()
    return metadata


_NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")
