from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import urlretrieve
import zipfile


LCB25_BASE_URL = "https://www.skies-village.it/webtools/Databases/LCB25/"
LCB25_DATASETS = ("PCS2", "SE", "HPCS2")


@dataclass(frozen=True)
class LCB25Dataset:
    label: str
    url: str
    archive_name: str


def lcb25_dataset_url(label: str) -> str:
    normalized = str(label).strip().upper()
    if normalized not in LCB25_DATASETS:
        raise ValueError(f"unsupported LCB25 dataset {label!r}; expected one of {LCB25_DATASETS}")
    return urljoin(LCB25_BASE_URL, f"{normalized}.zip")


def lcb25_download_plan() -> tuple[LCB25Dataset, ...]:
    return tuple(
        LCB25Dataset(label=label, url=lcb25_dataset_url(label), archive_name=f"{label}.zip")
        for label in LCB25_DATASETS
    )


def download_lcb25_dataset(label: str, target_dir: Path) -> Path:
    """Download one LCB25 ZIP archive to `target_dir`.

    Network use is explicit; tests cover URL planning without downloading.
    """
    dataset = LCB25Dataset(
        label=str(label).strip().upper(),
        url=lcb25_dataset_url(label),
        archive_name=f"{str(label).strip().upper()}.zip",
    )
    outdir = Path(target_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    archive = outdir / dataset.archive_name
    urlretrieve(dataset.url, archive)
    return archive


def extract_lcb25_archive(archive: Path, target_dir: Path) -> tuple[Path, ...]:
    """Extract an LCB25 archive and return extracted XYZ files."""
    archive = Path(archive)
    outdir = Path(target_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as zf:
        for member in zf.infolist():
            destination = (outdir / member.filename).resolve()
            if outdir.resolve() not in destination.parents and destination != outdir.resolve():
                raise ValueError(f"unsafe path in LCB25 archive: {member.filename}")
        zf.extractall(outdir)
    return tuple(sorted(path for path in outdir.rglob("*.xyz") if path.is_file()))
