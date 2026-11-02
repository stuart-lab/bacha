"""Bacha — bias-adjusted chromatin accessibility scoring for single-cell ATAC-seq."""

from importlib.metadata import PackageNotFoundError, version

from ._count import count
from ._summary import summary

try:
    __version__ = version("bacha")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = ["count", "summary", "__version__"]
