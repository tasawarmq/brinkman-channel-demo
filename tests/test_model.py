import unittest

import numpy as np

from brinkman_channel import flow_rate, velocity


PARAMS = {"half_width": 1.3, "pressure_gradient": 2.0, "viscosity": 0.8}


class ChannelFlowTests(unittest.TestCase):
    def test_no_slip_and_symmetry(self):
        h = PARAMS["half_width"]
        y = np.array([-h, -h / 3, 0.0, h / 3, h])
        for alpha in [0.0, 1e-8, 0.1, 1.0, 3.0, 10.0, 1000.0]:
            with self.subTest(alpha=alpha):
                u = velocity(y, alpha=alpha, **PARAMS)
                np.testing.assert_allclose(u[[0, -1]], 0.0, atol=1e-12)
                np.testing.assert_allclose(u, u[::-1], rtol=0, atol=1e-13)

    def test_zero_drag_limit(self):
        h, g, mu = PARAMS["half_width"], PARAMS["pressure_gradient"], PARAMS["viscosity"]
        y = np.linspace(-h, h, 11)
        expected = g * (h**2 - y**2) / (2 * mu)
        np.testing.assert_allclose(velocity(y, alpha=0.0, **PARAMS), expected, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(velocity(y, alpha=1e-7, **PARAMS), expected, rtol=1e-12, atol=1e-12)
        self.assertAlmostEqual(flow_rate(alpha=0.0, **PARAMS), 2 * g * h**3 / (3 * mu))

    def test_analytic_flux_matches_quadrature(self):
        h = PARAMS["half_width"]
        y = np.linspace(-h, h, 20001)
        for alpha in [0.0, 0.001, 0.1, 1.0, 3.0, 8.0]:
            with self.subTest(alpha=alpha):
                numerical_flux = np.trapezoid(velocity(y, alpha=alpha, **PARAMS), y)
                analytic_flux = flow_rate(alpha=alpha, **PARAMS)
                self.assertLessEqual(abs(numerical_flux / analytic_flux - 1), 2e-7)

    def test_flux_decreases_with_screening(self):
        rates = [flow_rate(alpha=alpha, **PARAMS) for alpha in [0, 0.01, 0.1, 0.5, 1, 3, 10]]
        self.assertTrue(all(left > right for left, right in zip(rates, rates[1:])))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            flow_rate(half_width=-1, pressure_gradient=1, viscosity=1, alpha=1)
        with self.assertRaises(ValueError):
            velocity(np.array([1.01]), half_width=1, pressure_gradient=1, viscosity=1, alpha=1)


if __name__ == "__main__":
    unittest.main()
