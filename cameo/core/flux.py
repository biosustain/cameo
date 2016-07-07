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

from __future__ import absolute_import, print_function

from math import floor, ceil, trunc
from numbers import Real, Number


class Flux(Real):
    def __init__(self, value, precision=1e-6):
        self.value = value
        self._precision = precision
        self._ndecimals = None

    @property
    def precision(self):
        return self._precision

    @precision.setter
    def precision(self, precision):
        self._precision = 0
        self._ndecimals = None

    @property
    def ndecimals(self):
        if self._ndecimals is None:
            self._ndecimals = self.__ndecimals__()
        return self._ndecimals

    def __ndecimals__(self):
        number = 0.1
        if self._precision > number:
            return 0
        else:
            i = 1
            while self._precision < number:
                number /= 10
                i += 1

            return i

    def __round__(self, n=None):
        return round(self.value, n)

    def __truediv__(self, other):
        if isinstance(other, Flux):
            value = self.value / other.real
            precision = max(self.precision, other.precision)
        elif isinstance(other, (Number, float, int)):
            value = self.value / other
            precision = self.precision
        else:
            raise ValueError(other)

        return Flux(value, precision=precision)

    def __add__(self, other):
        if isinstance(other, Flux):
            value = self.value + other.real
            precision = max(self.precision, other.precision)
        elif isinstance(other, (Number, float, int)):
            value = self.value + other
            precision = self.precision
        else:
            raise ValueError(other)

        return Flux(value, precision=precision)

    def __iadd__(self, other):
        if isinstance(other, Flux):
            value = self.value + other.real
            precision = max(self.precision, other.precision)
        elif isinstance(other, (Number, float, int)):
            value = self.value + other
            precision = self.precision
        else:
            raise ValueError(other)

        self.value = value
        self.precision = precision

    def __mul__(self, other):
        if isinstance(other, Flux):
            value = self.value * other.real
            precision = max(self.precision, other.precision)
        elif isinstance(other, (Number, float, int)):
            value = self.value * other
            precision = self.precision
        else:
            raise ValueError(other)

        return Flux(value, precision)

    def __imul__(self, other):
        if isinstance(other, Flux):
            value = self.value * other.real
            precision = max(self.precision, other.precision)
        elif isinstance(other, (Number, float, int)):
            value = self.value * other
            precision = self.precision
        else:
            raise ValueError(other)

        self.value = value
        self.precision = precision

    def __pow__(self, exponent):
        return Flux(self.value**exponent, self.precision)

    def __rtruediv__(self, other):
        return self / other

    def __abs__(self):
        return Flux(abs(self.value), self.precision)

    def __neg__(self):
        return Flux(-self.value, self.precision)

    @property
    def real(self):
        return round(self.value, self.ndecimals)

    def __complex__(self):
        return complex(self.real)

    def __rpow__(self, base):
        return base ** self.real

    def __rmul__(self, other):
        return self * other

    def __pos__(self):
        return Flux(self.value, self.precision)

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return -self + other

    def __str__(self):
        return "%f +- %f" % (self.real, self.precision)

    def __eq__(self, other):
        if isinstance(other, Flux):
            value = self.value - other.real
            precision = max(self.precision, other.precision)
        elif isinstance(other, (Number, float, int)):
            value = self.value - other
            precision = self.precision
        else:
            return False
        return abs(value) <= precision

    def __gt__(self, other):
        if isinstance(other, Flux):
            value = self.value - other.real
            precision = max(self.precision, other.precision)
        elif isinstance(other, (Number, float, int)):
            value = self.value - other
            precision = self.precision
        else:
            return False

        return (value + precision) > 0 and (value - precision) > 0

    def __lt__(self, other):
        if isinstance(other, Flux):
            value = self.value - other.real
            precision = max(self.precision, other.precision)
        elif isinstance(other, (Number, float, int)):
            value = self.value - other
            precision = self.precision
        else:
            return False

        return (value + precision) < 0

    def __ge__(self, other):
        return (self > other) or (self == other)

    def __le__(self, other):
        return (self < other) or (self == other)

    def __ne__(self, other):
        if isinstance(other, Flux):
            value = self.value - other.real
            precision = max(self.precision, other.precision)
        elif isinstance(other, (Number, float, int)):
            value = self.value - other
            precision = self.precision
        else:
            return False
        return abs(value) > precision

    def __repr__(self):
        return str(self)

    def __mod__(self, other):
        return self.real % other

    def __floordiv__(self, other):
        floor(self / other)

    def __rfloordiv__(self, other):
        return floor(other / self)

    def __trunc__(self):
        return trunc(self.real)

    def __floor__(self):
        return floor(self.value + self.precision)

    def __float__(self):
        return self.real

    def __rmod__(self, other):
        return other % self.real

    def __ceil__(self):
        return ceil(self.value - self.precision)
