"""Pressure-driven, plane Brinkman flow with no-slip walls.

The dimensional equation is
    mu*u''(y) - mu*alpha**2*u(y) = -pressure_gradient
on -half_width <= y <= half_width, with u(+-half_width) = 0.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _validate(half_width: float, viscosity: float, alpha: float) -> None:
    if not np.isfinite(half_width) or half_width <= 0:
        raise ValueError("half_width must be finite and positive")
    if not np.isfinite(viscosity) or viscosity <= 0:
        raise ValueError("viscosity must be finite and positive")
    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and nonnegative")


def velocity(
    y: ArrayLike,
    *,
    half_width: float,
    pressure_gradient: float,
    viscosity: float,
    alpha: float,
) -> NDArray[np.float64]:
    """Return axial velocity for y in [-half_width, half_width].

    `pressure_gradient` means -dp/dx, positive for flow in +x.
    The small-screening expansion avoids a 0/0 cancellation.
    """
    _validate(half_width, viscosity, alpha)
    if not np.isfinite(pressure_gradient):
        raise ValueError("pressure_gradient must be finite")
    y_array = np.asarray(y, dtype=float)
    if not np.all(np.isfinite(y_array)):
        raise ValueError("y must contain finite values")
    if np.any(np.abs(y_array) > half_width * (1 + 1e-12)):
        raise ValueError("y values must lie between the walls")

    x = y_array / half_width
    beta = alpha * half_width
    if not np.isfinite(beta):
        raise ValueError("alpha * half_width must be finite")
    scale = pressure_gradient * half_width**2 / viscosity

    if beta < 1e-3:
        result = (1 - x**2) / 2 + beta**2 * (-5 + 6 * x**2 - x**4) / 24
    else:
        # Avoid cosh overflow when the screening parameter is large.
        ratio = (np.exp(-beta * (1 - x)) + np.exp(-beta * (1 + x))) / (1 + np.exp(-2 * beta))
        result = (1 - ratio) / beta**2
    return np.asarray(scale * result, dtype=float)


def flow_rate(
    *,
    half_width: float,
    pressure_gradient: float,
    viscosity: float,
    alpha: float,
) -> float:
    """Return volumetric flux per unit depth perpendicular to the x-y plane."""
    _validate(half_width, viscosity, alpha)
    if not np.isfinite(pressure_gradient):
        raise ValueError("pressure_gradient must be finite")
    beta = alpha * half_width
    if not np.isfinite(beta):
        raise ValueError("alpha * half_width must be finite")
    scale = pressure_gradient * half_width**3 / viscosity
    if beta < 1e-2:
        dimensionless_flux = 2 / 3 - 4 * beta**2 / 15 + 34 * beta**4 / 315
    else:
        dimensionless_flux = 2 * (1 - np.tanh(beta) / beta) / beta**2
    return float(scale * dimensionless_flux)
