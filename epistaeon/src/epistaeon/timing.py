"""
epistaeon/timing.py
-------------------
Where along a branch each substitution falls.

Deliberately limited. The relative position of a substitution is reported from
the phi increments, which the data supports. Absolute times in generations
require an effective population size and a selection coefficient, which nothing
in this repository provides, so they are opt-in and guarded.

The white paper's sweep expression, t_sweep ~ 2*ln(2*Ne*s)/s, is undefined for
the near-neutral regime its own permissive-drift argument invokes: when
2*Ne*s <= 1 the logarithm is non-positive and the formula returns a negative
time. We raise instead of returning nonsense.
"""

from dataclasses import dataclass
from math import log
from typing import Dict, List, Optional, Sequence

import numpy as np


class TimingNotSupported(ValueError):
    """Raised when absolute timing is requested outside its valid regime."""


@dataclass
class RelativePosition:
    """Fractional position along the branch, from phi increments."""

    unit_names: List[str]
    order: List[int]
    phi_at_step: List[float]
    fraction: List[float]
    monotone: bool = True
    note: str = ""

    def as_dict(self) -> Dict[str, float]:
        return {
            self.unit_names[u]: float(self.fraction[k])
            for k, u in enumerate(self.order)
        }


def relative_positions(
    unit_names: Sequence[str],
    order: Sequence[int],
    phi_values: Sequence[float],
) -> RelativePosition:
    """Place each substitution by the share of total phi travelled before it.

    `phi_values` is phi at each state along the ordering, starting at the
    ancestor and ending at the derived sequence (len = len(order) + 1).
    """
    phi = np.asarray(phi_values, dtype=float)
    if len(phi) != len(order) + 1:
        raise ValueError("phi_values must have one more entry than order")
    span = float(phi[-1] - phi[0])
    if abs(span) < 1e-12:
        frac = [float(k + 1) / len(order) for k in range(len(order))]
    else:
        frac = [float((phi[k + 1] - phi[0]) / span) for k in range(len(order))]
    monotone = all(b >= a - 1e-9 for a, b in zip(frac, frac[1:]))
    note = "" if monotone else (
        "phi is not monotone along this ordering, so a fraction may fall "
        "outside [0, 1]. phi measures divergence from the anchor, which can "
        "overshoot the derived sequence; read these as phi levels, not as "
        "positions in time."
    )
    return RelativePosition(
        unit_names=list(unit_names), order=list(order),
        phi_at_step=[float(v) for v in phi], fraction=frac,
        monotone=monotone, note=note,
    )


def drift_time_generations(ne: float) -> float:
    """Expected fixation time of a neutral mutation, ~4*Ne generations."""
    if ne <= 0:
        raise TimingNotSupported(f"Ne must be positive (got {ne})")
    return 4.0 * float(ne)


def sweep_time_generations(ne: float, s: float) -> float:
    """Sweep duration ~2*ln(2*Ne*s)/s, guarded against its invalid regime."""
    if ne <= 0:
        raise TimingNotSupported(f"Ne must be positive (got {ne})")
    if s <= 0:
        raise TimingNotSupported(
            f"selection coefficient must be positive for a sweep (got {s})"
        )
    product = 2.0 * float(ne) * float(s)
    if product <= 1.0:
        raise TimingNotSupported(
            f"2*Ne*s = {product:.4g} <= 1: the mutation is effectively neutral, "
            "where the sweep expression is undefined and would return a negative "
            "time. Use the drift timescale instead."
        )
    return 2.0 * log(product) / float(s)


def absolute_timing(
    positions: RelativePosition,
    ne: Optional[float],
    s: Optional[float],
    generation_time_years: Optional[float] = None,
) -> Dict[str, Dict[str, float]]:
    """Convert relative positions to timescales. Requires Ne and s explicitly."""
    if ne is None or s is None:
        raise TimingNotSupported(
            "absolute timing needs both --ne and --s; no source in this "
            "repository provides them, so they must be supplied deliberately"
        )
    drift = drift_time_generations(ne)
    sweep = sweep_time_generations(ne, s)
    out: Dict[str, Dict[str, float]] = {}
    for name, frac in positions.as_dict().items():
        entry = {
            "fraction_of_branch": frac,
            "drift_timescale_generations": drift,
            "sweep_timescale_generations": sweep,
        }
        if generation_time_years:
            entry["drift_timescale_years"] = drift * float(generation_time_years)
            entry["sweep_timescale_years"] = sweep * float(generation_time_years)
        out[name] = entry
    return out
