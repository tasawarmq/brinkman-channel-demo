# Brinkman flow between parallel plates

Python programs for pressure-driven Brinkman flow between two parallel no-slip plates.

The analytical solution is used to calculate velocity profiles and flow rates for different resistance parameters.

## Running the code

Install:

```bash
python -m pip install -e .
```

Run the tests:

```bash
python -m unittest discover -s tests -v
```

Calculate flow rates:

```bash
brinkman-sweep --alphas 0,0.1,0.5,1,2,5 --output results/flow_rates.csv
```

MIT License.
