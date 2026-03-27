"""Custom exceptions for PhytClust."""


class PhytClustError(Exception):
    """Base exception for all PhytClust errors."""
    pass


class ValidationError(PhytClustError):
    """Raised when input validation fails."""
    pass


class ConfigurationError(PhytClustError):
    """Raised when configuration parameters are invalid."""
    pass


class ComputationError(PhytClustError):
    """Raised when an internal computation fails."""
    pass


class DataError(PhytClustError):
    """Raised when data is malformed or missing."""
    pass


class InvalidKError(ConfigurationError):
    """Raised when k value is invalid."""
    pass


class InvalidTreeError(ValidationError):
    """Raised when tree structure is invalid."""
    pass


class MissingDPTableError(ComputationError):
    """Raised when DP table is not computed but required."""
    pass


class InvalidClusteringError(ComputationError):
    """Raised when clustering cannot be performed."""
    pass
