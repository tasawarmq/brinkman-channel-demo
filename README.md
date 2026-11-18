# Brinkman flow between parallel plates

A small, self-contained Python example of steady pressure-driven flow between two no-slip plates with linear Brinkman drag.

This is an **independently written educational demonstration**, not the source code or a reproduction of any published manuscript. It contains no unpublished results, collaborator files, experimental data or university-specific software.

## Model

For a fluid between plates at `y = -H` and `y = H`, take the pressure-gradient magnitude `G = -dp/dx > 0`, dynamic viscosity `mu`, and inverse screening length `alpha >= 0`. The axial velocity satisfies

```text
mu * d²u/dy² - mu * alpha² * u = -G,
u(-H) = u(H) = 0.
```

The solution for `alpha > 0` is

```text
u(y) = G / (mu * alpha²) * [1 - cosh(alpha*y)/cosh(alpha*H)].
```

For `alpha = 0` it reduces to plane Poiseuille flow,

```text
u(y) = G / (2*mu) * (H²-y²).
```

The volumetric flow rate per unit span is

```text
Q = 2*G/(mu*alpha²) * [H - tanh(alpha*H)/alpha]
```

for `alpha > 0`, with `Q = 2*G*H³/(3*mu)` at `alpha = 0`.

Small-`alpha` series expansions are used to avoid cancellation when evaluating the formulas numerically.

## Installation and use

Requires Python 3.10+ and NumPy. From this directory:

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
brinkman-sweep --alphas 0,0.1,0.5,1,2,5 --output results/flow_rates.csv
```

The output CSV contains an `alpha` value, the dimensionless screening parameter `alpha_H`, the computed flux, and the flux divided by its zero-drag value.

## Reproducibility checks

- The velocity is zero on each plate and symmetric about the midplane.
- The small-screening limit agrees with Poiseuille flow.
- Numerical integration of the velocity agrees with the analytical flux.
- At fixed pressure gradient, increasing drag reduces the volume flux.

## Scope and attribution

This is a standard textbook Brinkman channel problem, not a new research result. It has been newly implemented as a compact sample of numerical modelling, parameter sweeps and verification. It does not represent the author's original choanoflagellate, flagellar motor, or droplet-interface code.

Released under the MIT License (see [LICENSE](LICENSE)).
