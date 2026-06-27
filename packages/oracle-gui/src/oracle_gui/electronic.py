from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from oracle_core import has_section, read_sectioned_lines, section_content

from .commands import (
    MORBVIS_URL,
    OracleGuiCommand,
    external_viewer_command,
    molden_command,
    morbvis_command,
)


ELECTRONIC_SECTION_NAMES = ("ELECTRONIC", "TRANSITIONS", "ORBITALS")


@dataclass(frozen=True)
class ElectronicTable:
    title: str
    columns: tuple[str, ...]
    rows: tuple[tuple[str, ...], ...]


@dataclass(frozen=True)
class ElectronicGuiState:
    xyzin: Path
    exists: bool
    ready: bool
    sections: ElectronicTable
    transitions: ElectronicTable
    orbitals: ElectronicTable
    messages: tuple[str, ...] = ()


class OracleElectronicController:
    def __init__(self, xyzin: Path | str | None = None) -> None:
        self.xyzin = None if xyzin is None else Path(xyzin)

    def set_xyzin(self, xyzin: Path | str | None) -> ElectronicGuiState | None:
        self.xyzin = None if xyzin is None else Path(xyzin)
        if self.xyzin is None:
            return None
        return self.state()

    def state(self) -> ElectronicGuiState:
        if self.xyzin is None:
            raise ValueError("no ORACLE xyzin project is loaded")
        return load_electronic_gui_state(self.xyzin)

    def molden_command(
        self,
        target: Path | str,
        *,
        executable: str = "molden",
    ) -> OracleGuiCommand:
        return molden_command(target, executable=executable)

    def avogadro_command(
        self,
        target: Path | str,
        *,
        executable: str = "avogadro2",
    ) -> OracleGuiCommand:
        return external_viewer_command(target, executable=executable, label="Open in Avogadro")

    def morbvis_command(self, *, url: str = MORBVIS_URL) -> OracleGuiCommand:
        return morbvis_command(url=url)


def load_electronic_gui_state(path: Path | str) -> ElectronicGuiState:
    target = Path(path)
    if not target.exists():
        return ElectronicGuiState(
            xyzin=target,
            exists=False,
            ready=False,
            sections=ElectronicTable("Sections", ("Section", "Present", "Lines"), ()),
            transitions=ElectronicTable("Transitions", ("Line", "Record"), ()),
            orbitals=ElectronicTable("Orbitals", ("Line", "Record"), ()),
            messages=(f"Missing file: {target}",),
        )
    section_rows = []
    for name in ELECTRONIC_SECTION_NAMES:
        content = section_content(read_sectioned_lines(target), name) if has_section(target, name) else []
        section_rows.append((name, "yes" if content else "no", str(len(content))))
    ready = any(row[1] == "yes" for row in section_rows)
    messages = () if ready else ("missing #ELECTRONIC/#TRANSITIONS/#ORBITALS sections",)
    return ElectronicGuiState(
        xyzin=target,
        exists=True,
        ready=ready,
        sections=ElectronicTable("Sections", ("Section", "Present", "Lines"), tuple(section_rows)),
        transitions=_content_table(target, "TRANSITIONS"),
        orbitals=_content_table(target, "ORBITALS"),
        messages=messages,
    )


def electronic_gui_state_lines(state: ElectronicGuiState) -> list[str]:
    return [
        f"xyzin: {state.xyzin}",
        f"ready: {int(state.ready)}",
        *state.messages,
    ]


def _content_table(target: Path, section: str) -> ElectronicTable:
    if not has_section(target, section):
        return ElectronicTable(section.title(), ("Line", "Record"), ())
    content = [
        line.strip()
        for line in section_content(read_sectioned_lines(target), section)
        if line.strip() and not line.strip().upper().startswith("SCHEMA ")
    ]
    parsed = _columns_table(section, content)
    if parsed is not None:
        return parsed
    return ElectronicTable(
        section.title(),
        ("Line", "Record"),
        tuple((str(index), line) for index, line in enumerate(content, start=1)),
    )


def _columns_table(section: str, content: list[str]) -> ElectronicTable | None:
    columns: tuple[str, ...] | None = None
    rows: list[tuple[str, ...]] = []
    for raw in content:
        if raw.upper().startswith("COLUMNS "):
            columns = tuple(raw.split()[1:])
            continue
        if columns is None:
            continue
        parts = tuple(raw.split())
        if len(parts) == len(columns):
            rows.append(parts)
    if columns is None:
        return None
    return ElectronicTable(section.title(), columns, tuple(rows))
