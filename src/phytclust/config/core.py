from dataclasses import dataclass
from typing import Literal


@dataclass
class CoreConfig:
    """Core algorithm behavior toggles for PhytClust."""

    optimize_polytomies: bool = True
    polytomy_mode: Literal["hard", "soft"] = "hard"
    soft_polytomy_max_degree: int = 18
    no_split_zero_length: bool = False
    preserve_dp_tables: bool = False
