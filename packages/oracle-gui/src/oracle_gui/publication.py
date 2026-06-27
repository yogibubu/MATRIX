from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import io

from oracle_thermo import read_thermo_section
from oracle_thermo.models import THERMO_KEYS, THERMO_LABELS, ThermoSection


@dataclass(frozen=True)
class PublicationExportResult:
    outdir: Path
    paths: tuple[Path, ...]


def export_thermo_table(
    xyzin: Path | str,
    outdir: Path | str,
    *,
    basename: str | None = None,
    formats: tuple[str, ...] = ("csv", "tex", "svg", "pdf"),
) -> PublicationExportResult:
    target = Path(xyzin)
    section = read_thermo_section(target)
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    stem = basename or f"{target.stem}.thermo"
    requested = tuple(_normalize_format(item) for item in formats)
    rows = _thermo_rows(section)
    paths: list[Path] = []
    for fmt in requested:
        if fmt == "csv":
            paths.append(_write_csv(destination / f"{stem}.csv", rows))
        elif fmt == "tex":
            paths.append(_write_latex(destination / f"{stem}.tex", rows))
        elif fmt == "svg":
            paths.append(_write_svg(destination / f"{stem}.svg", rows))
        elif fmt == "pdf":
            paths.append(_write_pdf(destination / f"{stem}.pdf", rows))
        else:
            raise ValueError(f"unsupported publication export format: {fmt}")
    return PublicationExportResult(destination, tuple(paths))


def _normalize_format(value: str) -> str:
    text = value.strip().lower().lstrip(".")
    if text == "latex":
        return "tex"
    return text


def _thermo_rows(section: ThermoSection) -> list[list[str]]:
    rows = [["component", *THERMO_KEYS]]
    for label in THERMO_LABELS:
        contribution = section.contribution(label)
        if contribution is None:
            continue
        rows.append(
            [
                label,
                *[_format_value(getattr(contribution, key)) for key in THERMO_KEYS],
            ]
        )
    return rows


def _format_value(value: float | None) -> str:
    return "" if value is None else f"{float(value):.12g}"


def _write_csv(path: Path, rows: list[list[str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)
    return path


def _write_latex(path: Path, rows: list[list[str]]) -> Path:
    column_spec = "l" + "r" * (len(rows[0]) - 1)
    lines = [
        rf"\begin{{tabular}}{{{column_spec}}}",
        r"\hline",
        " & ".join(_latex_cell(cell) for cell in rows[0]) + r" \\",
        r"\hline",
    ]
    for row in rows[1:]:
        lines.append(" & ".join(_latex_cell(cell) for cell in row) + r" \\")
    lines.extend([r"\hline", r"\end{tabular}"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_svg(path: Path, rows: list[list[str]]) -> Path:
    width = 1180
    line_height = 24
    height = 40 + line_height * len(rows)
    columns = len(rows[0])
    col_width = width / columns
    text = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Helvetica,Arial,sans-serif;font-size:13px}'
        '.head{font-weight:bold}</style>',
    ]
    for row_index, row in enumerate(rows):
        y = 28 + row_index * line_height
        cls = ' class="head"' if row_index == 0 else ""
        for column_index, cell in enumerate(row):
            x = 12 + column_index * col_width
            text.append(f'<text{cls} x="{x:.1f}" y="{y:.1f}">{_xml_escape(cell)}</text>')
    text.append("</svg>")
    path.write_text("\n".join(text) + "\n", encoding="utf-8")
    return path


def _write_pdf(path: Path, rows: list[list[str]]) -> Path:
    lines = ["  ".join(row) for row in rows]
    stream_lines = ["BT", "/F1 8 Tf", "40 780 Td"]
    for line in lines:
        stream_lines.append(f"({_pdf_escape(line)}) Tj")
        stream_lines.append("0 -14 Td")
    stream_lines.append("ET")
    stream = "\n".join(stream_lines).encode("latin-1", errors="replace")
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] "
        b"/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        b"<< /Length " + str(len(stream)).encode("ascii") + b" >>\nstream\n" + stream + b"\nendstream",
    ]
    buffer = io.BytesIO()
    buffer.write(b"%PDF-1.4\n")
    offsets = [0]
    for index, obj in enumerate(objects, start=1):
        offsets.append(buffer.tell())
        buffer.write(f"{index} 0 obj\n".encode("ascii"))
        buffer.write(obj)
        buffer.write(b"\nendobj\n")
    xref = buffer.tell()
    buffer.write(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    buffer.write(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        buffer.write(f"{offset:010d} 00000 n \n".encode("ascii"))
    buffer.write(
        f"trailer << /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n".encode(
            "ascii"
        )
    )
    path.write_bytes(buffer.getvalue())
    return path


def _latex_cell(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
    }
    return "".join(replacements.get(char, char) for char in text)


def _xml_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _pdf_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
