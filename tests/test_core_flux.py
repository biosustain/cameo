# Copyright 2016 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import unittest

from cameo.core.flux import Flux
from math import ceil, floor


class FluxTestCase(unittest.TestCase):
    def setUp(self):
        self.fluxes = [
            Flux(10),
            Flux(5),
            Flux(1e-5),
            Flux(1e-6),
            Flux(1e-7),
            Flux(0),
            Flux(-1e-7),
            Flux(-1e-6),
            Flux(-1e-5),
            Flux(-5),
            Flux(-10)
        ]

    def test_add(self):
        self.assertEqual(self.fluxes[0] + self.fluxes[1], Flux(15))
        self.assertEqual(self.fluxes[0] + 5, Flux(15))
        self.assertEqual(10 + self.fluxes[1], Flux(15))

        less_precise_flux = Flux(10, precision=1e-5)
        flux_sum = self.fluxes[0] + less_precise_flux
        self.assertEqual(flux_sum.precision, less_precise_flux.precision)

    def test_ndecimals(self):
        self.assertTrue(self.fluxes[0].ndecimals, 6)

        less_precise_flux = Flux(10, precision=1e-5)
        more_precise_flux = Flux(10, precision=1e-10)

        self.assertTrue(less_precise_flux.ndecimals, 5)
        self.assertTrue(more_precise_flux.ndecimals, 5)

    def test_sub(self):
        self.assertEqual(self.fluxes[0] - self.fluxes[1], Flux(5))
        self.assertEqual(self.fluxes[0] - 5, Flux(5))
        self.assertEqual(10 - self.fluxes[1], Flux(5))

        less_precise_flux = Flux(10, precision=1e-5)
        flux_sum = self.fluxes[0] - less_precise_flux
        self.assertEqual(flux_sum.precision, less_precise_flux.precision)

    def test_multiply(self):
        self.assertEqual(self.fluxes[0] * self.fluxes[1], 50)
        self.assertEqual(self.fluxes[1] * self.fluxes[0], 50)
        self.assertEqual(self.fluxes[1] * self.fluxes[0], Flux(50))
        self.assertEqual(self.fluxes[0] * self.fluxes[5], 0)

        less_precise_flux = Flux(10, precision=1e-5)
        self.assertEqual((less_precise_flux * self.fluxes[5]).precision, 1e-5)

    def test_truediv(self):
        self.assertEqual(self.fluxes[0] / self.fluxes[1], 2)
        self.assertEqual(self.fluxes[1] / self.fluxes[0], 0.5)
        self.assertEqual(self.fluxes[1] / self.fluxes[0], Flux(0.5))
        with self.assertRaises(ZeroDivisionError):
            self.fluxes[0].__truediv__(self.fluxes[5])

        less_precise_flux = Flux(10, precision=1e-5)
        self.assertEqual((less_precise_flux / self.fluxes[1]).precision, 1e-5)

    def test_exponential(self):
        self.assertEqual(self.fluxes[0] ** 2, 10 ** 2)
        self.assertEqual(self.fluxes[1] ** 2, 5 ** 2)
        self.assertEqual(self.fluxes[3] ** 2, 1e-5 ** 2)
        self.assertEqual(Flux(0.0000011, precision=1e-5) ** 2, 1e-6 ** 2)

    def test_comparator(self):
        more_precise_flux = Flux(10, precision=1e-10)

        self.assertEqual(self.fluxes[0], more_precise_flux)
        self.assertEqual(self.fluxes[0], 10+1e-6)
        self.assertEqual(self.fluxes[0], 10-1e-6)
        self.assertNotEqual(self.fluxes[0], 10+1e-5)
        self.assertNotEqual(self.fluxes[0], 10-1e-5)
        self.assertGreater(self.fluxes[0], 9)
        self.assertGreaterEqual(self.fluxes[0], 9)
        self.assertEqual(self.fluxes[0], 10)
        self.assertGreaterEqual(self.fluxes[0], 10)

        self.assertLess(self.fluxes[0], 11)
        self.assertLessEqual(self.fluxes[0], 11)
        self.assertLessEqual(self.fluxes[0], 10)

    def test_bool(self):
        self.assertTrue(self.fluxes[0])
        self.assertTrue(self.fluxes[1])
        self.assertTrue(self.fluxes[2])
        self.assertFalse(self.fluxes[3])
        self.assertFalse(self.fluxes[4])
        self.assertFalse(self.fluxes[5])
        self.assertFalse(self.fluxes[6])
        self.assertFalse(self.fluxes[7])
        self.assertTrue(self.fluxes[8])
        self.assertTrue(self.fluxes[9])
        self.assertTrue(self.fluxes[10])

    def test_ceil_floor(self):
        self.assertEqual(ceil(self.fluxes[0]), 10)
        self.assertEqual(floor(self.fluxes[0]), 10)
        self.assertEqual(ceil(self.fluxes[0] + 1e-6), 10)
        self.assertEqual(ceil(self.fluxes[0] + 1e-5), 11)
        self.assertEqual(floor(self.fluxes[0] + 1e-6), 10)
        self.assertEqual(floor(self.fluxes[0] + 1e-5), 10)

