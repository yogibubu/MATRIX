from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import os
from pathlib import Path


GIC_CORPUS_ENV = "ORACLE_GIC_CORPUS"

ROLE_BY_SUFFIX = {
    ".inp": "legacy_gic_input",
    ".gjf": "gaussian_input",
    ".gau": "gaussian_input",
    ".fchk": "gaussian_fchk",
    ".log": "qm_output",
    ".out": "qm_output",
    ".msr": "morpheus_state",
    ".sum": "morpheus_summary",
    ".opt": "morpheus_optimized",
    ".form": "morpheus_formatted",
    ".vlt": "legacy_vlt",
    ".gbs": "basis_set",
    "": "unclassified",
}


@dataclass(frozen=True)
class GICCorpusEntry:
    path: Path
    name: str
    stem: str
    suffix: str
    role: str

    def record(self, *, root: Path | None = None) -> dict[str, str]:
        resolved_root = root.resolve() if root is not None else None
        resolved_path = self.path.resolve()
        if resolved_root is not None:
            try:
                display_path = str(resolved_path.relative_to(resolved_root))
            except ValueError:
                display_path = str(resolved_path)
        else:
            display_path = str(self.path)
        return {
            "name": self.name,
            "path": display_path,
            "stem": self.stem,
            "suffix": self.suffix,
            "role": self.role,
        }


@dataclass(frozen=True)
class GICCorpusSummary:
    root: Path
    entries: tuple[GICCorpusEntry, ...]

    @property
    def total_files(self) -> int:
        return len(self.entries)

    @property
    def suffix_counts(self) -> dict[str, int]:
        return dict(Counter(entry.suffix for entry in self.entries))

    @property
    def role_counts(self) -> dict[str, int]:
        return dict(Counter(entry.role for entry in self.entries))


class GICCorpusError(ValueError):
    """Raised when the GIC regression corpus cannot be discovered."""


def default_gic_corpus_root(repo_root: Path | None = None) -> Path:
    env_root = os.environ.get(GIC_CORPUS_ENV)
    if env_root:
        return Path(env_root).expanduser().resolve()
    if repo_root is None:
        repo_root = Path.cwd()
    return Path(repo_root).resolve() / "tests" / "fixtures" / "test_molecules" / "molecules"


def discover_gic_corpus(
    root: Path,
    *,
    suffixes: tuple[str, ...] | list[str] | None = None,
) -> tuple[GICCorpusEntry, ...]:
    target = Path(root).expanduser().resolve()
    if not target.is_dir():
        raise GICCorpusError(f"GIC regression corpus directory not found: {target}")

    requested = _normalize_suffixes(suffixes)
    entries: list[GICCorpusEntry] = []
    for path in sorted(item for item in target.iterdir() if item.is_file()):
        suffix = path.suffix.lower()
        if requested and suffix not in requested:
            continue
        entries.append(
            GICCorpusEntry(
                path=path,
                name=path.name,
                stem=path.stem,
                suffix=suffix,
                role=ROLE_BY_SUFFIX.get(suffix, "unclassified"),
            )
        )
    return tuple(entries)


def summarize_gic_corpus(
    root: Path,
    *,
    suffixes: tuple[str, ...] | list[str] | None = None,
) -> GICCorpusSummary:
    target = Path(root).expanduser().resolve()
    return GICCorpusSummary(root=target, entries=discover_gic_corpus(target, suffixes=suffixes))


def gic_corpus_records(
    summary: GICCorpusSummary,
    *,
    limit: int | None = None,
) -> list[dict[str, str]]:
    entries = _limited(summary.entries, limit)
    return [entry.record(root=summary.root) for entry in entries]


def format_gic_corpus_summary(summary: GICCorpusSummary) -> list[str]:
    lines = [
        f"ROOT {summary.root}",
        f"TOTAL_FILES {summary.total_files}",
    ]
    for suffix, count in sorted(summary.suffix_counts.items()):
        label = suffix or "<none>"
        lines.append(f"SUFFIX {label} {count}")
    for role, count in sorted(summary.role_counts.items()):
        lines.append(f"ROLE {role} {count}")
    return lines


def format_gic_corpus_paths(
    summary: GICCorpusSummary,
    *,
    limit: int | None = None,
) -> list[str]:
    return [str(entry.path) for entry in _limited(summary.entries, limit)]


def _limited(entries: tuple[GICCorpusEntry, ...], limit: int | None) -> tuple[GICCorpusEntry, ...]:
    if limit is None:
        return entries
    if limit < 0:
        raise GICCorpusError("corpus limit cannot be negative")
    return entries[:limit]


def _normalize_suffixes(suffixes: tuple[str, ...] | list[str] | None) -> set[str]:
    if not suffixes:
        return set()
    normalized: set[str] = set()
    for suffix in suffixes:
        text = suffix.strip().lower()
        if text == "<none>":
            text = ""
        elif text and not text.startswith("."):
            text = f".{text}"
        normalized.add(text)
    return normalized
