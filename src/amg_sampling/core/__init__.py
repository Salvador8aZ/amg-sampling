"""Mathematical core: pure Python, no third-party dependencies."""

from amg_sampling.core.configuration import TELOMERE, InitialConfiguration, Side, Telomere
from amg_sampling.core.matching import InvalidMatchingError, RejoinMatching

__all__ = [
    "TELOMERE",
    "InitialConfiguration",
    "InvalidMatchingError",
    "RejoinMatching",
    "Side",
    "Telomere",
]
