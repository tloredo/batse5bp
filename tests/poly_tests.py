import nose
from nose.tools import assert_almost_equal

from numpy import array, asarray, single, csingle, ones_like
from numpy.testing import (TestCase, assert_, assert_equal, assert_raises,
                           assert_array_equal, # assert_almost_equal,
                           run_module_suite)

from batse5bp.poly import *


# Collections of cases (defining parameters for a test & the expected result),
# checkers (for testing a single case), and generators applying checks to cases:

def check_poly_case(u, b, result):
    r = factor_times_poly(u, b)
    assert_array_equal(r, result)

poly_cases = [ # (u, b, result)
    (0, [1], [0, 1]),  # x*1 = x
    (1, [1], [-1, 1]),  # (x-1)*1 = -1 + x
    (1, [0, 1], [0, -1, 1]),  # (x-1)*x = 0 - x + x**2
    (3, [2, -3, 4], [-6, 11, -15, 4]),  # (x-3)*(4*x**2 -3*x + 2) = -6 + 11*x - 15*x**2 + 4*x**3
    ]

def test_poly_cases():
    for case in poly_cases:
        yield check_poly_case, case[0], array(case[1], dtype=float), \
                array(case[2], dtype=float)


def check_roots_case(uvals, result):
    r = roots2coefs(uvals)
    assert_array_equal(r, result)

roots_cases = [ # (uvals, result)
    ([0], [0, 1]),  # x-0 = x
    ([1], [-1, 1]),  # -1 + x
    ([1, -1], [-1, 0, 1]),  # (x-1)*(x+1) = -1 + x**2
    ([1, -1, -1, 1], [1, 0, -2, 0, 1]),  # (x**2-1)**2 = 1 - 2*x**2 + x**4
    ([2, -3], [-6, 1, 1]),  # (x-2)*(x+3) = -6 + x + x**2
    ([2, -3, 5], [30, -11, -4, 1]),  # (x-2)*(x+3)*(x-5) = 30 - 11*x - 4*x**2 + x**3
    ]

def test_roots_cases():
    for case in roots_cases:
        yield check_roots_case, case[0], array(case[1], dtype=float)


def check_intgl_case(coef, a, b, result):
    r = integrate_poly(coef, a, b)
    assert_almost_equal(r, result)

intgl_cases = [ # (coef, a, b, result)
    ([0], 0, 1, 0),  # \int dx 0 = 0
    ([1], 0, 1, 1),  # \int dx 1 = x
    ([1], 2, 5, 3.),  # \int dx 1 = x
    ([0, 2], -1, 1, 0),  # \int dx 2*x = x**2
    ([0, 2], 1, 3, 8),  # \int dx 2*x = x**2
    ([0, 2], -1, 3, 8),  # \int dx 2*x = x**2
    ([3, 4, 6], 0, 1, 7),  # \int dx 6*x**2 + 4*x + 3 = 2*x**3 + 2*x**2 + 3*x
    ([3, 4, 6], 2, 5, 285),  # \int dx 6*x**2 + 4*x + 3 = 2*x**3 + 2*x**2 + 3*x -> 315 - 30
    ]

def test_intgl_cases():
    for case in intgl_cases:
        yield check_intgl_case, array(case[0], dtype=float), case[1], case[2], case[3]
