"""TRINITY external-gradient geometry optimization contracts."""

from .workflow import (
    ORACLE_XYZ_TRINITY_SCHEMA,
    TRINITY_LONG_NAME,
    TrinitySection,
    parse_trinity_section,
    prepare_trinity_section,
    read_trinity_section,
    trinity_section_from_request,
    trinity_section_lines,
    trinity_section_summary_lines,
    write_trinity_manifest,
    write_trinity_section,
)

__all__ = [
    "ORACLE_XYZ_TRINITY_SCHEMA",
    "TRINITY_LONG_NAME",
    "TrinitySection",
    "parse_trinity_section",
    "prepare_trinity_section",
    "read_trinity_section",
    "trinity_section_from_request",
    "trinity_section_lines",
    "trinity_section_summary_lines",
    "write_trinity_manifest",
    "write_trinity_section",
]
