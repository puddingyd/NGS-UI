"""Byte-for-byte reproducible ``.npz`` files (numpy's own writer stamps the time)."""
from __future__ import annotations

import io
import zipfile
from pathlib import Path

import numpy as np

_STAMP = (1980, 1, 1, 0, 0, 0)


def save(path: Path, **arrays) -> None:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name in sorted(arrays):
            buf = io.BytesIO()
            np.lib.format.write_array(buf, np.asarray(arrays[name]), allow_pickle=False)
            info = zipfile.ZipInfo(f"{name}.npy", date_time=_STAMP)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, buf.getvalue())


def load(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}
