from __future__ import annotations

import zipfile

import pytest

from oracle_babel import (
    LCB25_DATASETS,
    extract_lcb25_archive,
    lcb25_dataset_url,
    lcb25_download_plan,
)


def test_lcb25_download_plan_uses_official_zip_names():
    plan = lcb25_download_plan()

    assert tuple(item.label for item in plan) == LCB25_DATASETS
    assert lcb25_dataset_url("PCS2").endswith("/PCS2.zip")
    assert lcb25_dataset_url("SE").endswith("/SE.zip")
    assert lcb25_dataset_url("HPCS2").endswith("/HPCS2.zip")


def test_lcb25_archive_extraction_returns_xyz_files(tmp_path):
    archive = tmp_path / "SE.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("water.xyz", "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n")
        zf.writestr("notes.txt", "ignore")

    extracted = extract_lcb25_archive(archive, tmp_path / "out")

    assert len(extracted) == 1
    assert extracted[0].name == "water.xyz"


def test_lcb25_archive_extraction_rejects_path_traversal(tmp_path):
    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("../escape.xyz", "1\nbad\nH 0 0 0\n")

    with pytest.raises(ValueError, match="unsafe path"):
        extract_lcb25_archive(archive, tmp_path / "out")
