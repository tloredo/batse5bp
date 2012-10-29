import nose
from nose.tools import assert_almost_equal

from numpy import asarray, single, csingle, ones_like
from numpy.testing import (TestCase, assert_, assert_equal, assert_raises,
                           assert_array_equal, # assert_almost_equal,
                           run_module_suite)

from batse5bp.prodquad import ProdQuad11, CompositeQuad


# Borrowed from NumPy tests:

# numpy_assert_almost_equal = assert_almost_equal

# def assert_almost_equal(a, b, **kw):
#     print 'assert:', a, b
#     if asarray(a).dtype.type in (single, csingle):
#         decimal = 6
#     else:
#         decimal = 12
#     numpy_assert_almost_equal(a, b, decimal=decimal, **kw)


# Function pairs (f,g) with known exact integrals, specified as
# attributes of the g factor as (f, a, b, result) tuples.

def f_1(x):
    return ones_like(x)

def f_x(x):
    return x

def g_x(x):
    return x
g_x.cases = [(f_1, 0, 1, .5), (f_x, 0., 1., 1/3.)]

def g_xm1(x):
    return x-1.
g_xm1.cases = [(f_1, 0, 1, -.5), (f_x, 0., 1., -1/6.)]

def g_x2(x):
    return x**2
g_x2.cases = [(f_1, 0, 1, 1/3.), (f_x, 0., 1., .25)]


# Checker and generators for tests vs. exact results:

def check_fg_case(pq, f, g, result):
    q = pq.quad_fg(f, g)
    assert_almost_equal(q, result)

trap_cases = [g_x, g_xm1]

def test_trap_cases():
    for g in trap_cases:
        for f, a, b, result in g.cases:
            pq = ProdQuad11(a, b, a, b, a, b)  # product trapezoid rule
            yield check_fg_case, pq, f, g, result

GL_cases = [g_x, g_xm1]

def test_GL_cases():
    for g in GL_cases:
        for f, a, b, result in g.cases:
            pq = ProdQuad11(a, b, a, b)  # g nodes will be Legendre roots
            yield check_fg_case, pq, f, g, result

arbnode_cases = [g_x, g_xm1]

def test_arbnode_cases():
    for g in arbnode_cases:
        for f, a, b, result in g.cases:
            pq = ProdQuad11(a, b, 0.2, .6, 0.3, 0.7)  # arbitrary nodes
            yield check_fg_case, pq, f, g, result


# Test composite rules via 2-element cases built from the single-rule cases
# by cutting [a,b] in two.

def check_comp_fg_case(cq, g, result):
    q = cq.quad(g)
    # The next print stmt strangely causes the test to fail.
    #print 'CQ range:', cq.l, cq.u, [(q.l, q.u) for q in cq.rules]
    #print q, result
    assert_almost_equal(q, result)

def test_comp_trap_cases():
    for g in trap_cases:
        for f, a, b, result in g.cases:
            m = b / 3.
            pq1 = ProdQuad11(a, m, a, m, a, m, f)  # product trapezoid rule
            pq2 = ProdQuad11(m, b, m, b, m, b)  # this one without f for init
            cq = CompositeQuad(pq1.quad_object(), pq2.quad_object(f))
            yield check_comp_fg_case, cq, g, result


# TODO:  Test non-ufunc signatures.


# This is used only if the script is directly run, e.g., "python prodquad_tests.py".
if __name__ == "__main__":
    # run_module_suite()  # NumPy's alternative to nose.main()

    # nose.main appears not to respect this option; use "nosetests -s".
    #argv = ['--nocapture']  # don't capture print statements--disable stdout capture plugin
    #nose.main(argv=argv)

    from numpy import *

    print
    print '*** Main***'

    a, b = 0, 1
    pq = ProdQuad11(a, b, a, b, a, b, f_x)
    print

    g = trap_cases[0]
    for f, a, b, result in g.cases:
        m = b / 2.
        pq1 = ProdQuad11(a, m, a, m, a, m, f)  # product trapezoid rule
        pq2 = ProdQuad11(m, b, m, b, m, b, f)  # this one without f for init
        q1, q2 = pq1.quad_object(), pq2.quad_object()
        cq = CompositeQuad(q1, q2)
        print f
        print cq.nodes
        print 'quads: ', q1.quad(g), q2.quad(g)
        print 'pquads:', pq1.quad_g(g), pq2.quad_fg(f, g)
        print cq.quad(g)
        print

