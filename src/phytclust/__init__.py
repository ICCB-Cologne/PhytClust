from .algo.core import PhytClust
from .exceptions import (
    PhytClustError,
    ValidationError,
    ConfigurationError,
    ComputationError,
    DataError,
    InvalidKError,
    InvalidTreeError,
    MissingDPTableError,
    InvalidClusteringError,
)
from importlib.metadata import version, PackageNotFoundError

__all__ = [
    "PhytClust",
    "PhytClustError",
    "ValidationError",
    "ConfigurationError",
    "ComputationError",
    "DataError",
    "InvalidKError",
    "InvalidTreeError",
    "MissingDPTableError",
    "InvalidClusteringError",
]

try:
    __version__ = version("phytclust")
except PackageNotFoundError:
    __version__ = "0.0.0"
