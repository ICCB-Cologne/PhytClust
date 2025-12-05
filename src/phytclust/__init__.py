from .algo.core import PhytClust
from importlib.metadata import version, PackageNotFoundError

__all__ = ["PhytClust"]

try:
    __version__ = version("phytclust")
except PackageNotFoundError:
    __version__ = "0.0.0"
