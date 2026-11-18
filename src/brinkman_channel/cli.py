"""Run a small parameter sweep and export a reproducible CSV table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .model import flow_rate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alphas", default="0,0.1,0.5,1,2,5")
    parser.add_argument("--half-width", type=float, default=1.0)
    parser.add_argument("--viscosity", type=float, default=1.0)
    parser.add_argument("--pressure-gradient", type=float, default=1.0)
    parser.add_argument("--output", type=Path, default=Path("results/flow_rates.csv"))
    args = parser.parse_args()

    alpha_values = [float(value.strip()) for value in args.alphas.split(",")]
    q0 = flow_rate(
        half_width=args.half_width,
        pressure_gradient=args.pressure_gradient,
        viscosity=args.viscosity,
        alpha=0,
    )
    if q0 == 0:
        parser.error("pressure gradient must be nonzero to compute relative flux")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["alpha", "alpha_H", "flux_per_unit_depth", "relative_flux"])
        for alpha in alpha_values:
            q = flow_rate(
                half_width=args.half_width,
                pressure_gradient=args.pressure_gradient,
                viscosity=args.viscosity,
                alpha=alpha,
            )
            writer.writerow([alpha, alpha * args.half_width, q, q / q0])

    print(f"Wrote {len(alpha_values)} rows to {args.output}")


if __name__ == "__main__":
    main()
