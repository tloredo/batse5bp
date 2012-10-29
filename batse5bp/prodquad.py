""" 
Basic product quadrature algorithms, for integrating the product of a
tabulated response function and a signal model.

These implement small symmetric and asymmetric interpolatory product
quadrature rules based on Lagrange polynomials, approximating

  mu = \int_a^b dx  f(x) g(x)

by

  \sum_{i=0}^m \sum_{j=0}^n  f(u_i) a_{ij} g(v_j)

In the intended applications, f(x) is a response function tabulated at a fixed
set of nodes (usually not equally spaced), and g(x) is a signal model that may
be freely evaluated anywhere.  Usually many integrals will be needed with the
same choice of f(x), but with various g(x) (i.e., signals with different
choices of model parameters).

{u_i} (m values) and {v_j} (n values) are nodes for the quadrature rule.

a_{ij} is a matrix of quadrature weights.

The full integral over [a,b] will be handled via a compound rule, applying
small-order (m,n) rules to sub-intervals determined by the tabulated
nodes for f(x).

Since f(x) is tabulated, the rules here have u_0=a and u_m=b.  For the m=1
case (two u_i nodes), the rules are thus trapezoidal for f(x).

For g(x), the nodes may be freely specified.  Using nodes located at
zeros of orthogonal polynomials leads to better convergence behavior
than regularly-spaced nodes.  Such nodes typically produce open rules
(i.e., with v_0 != a and v_n != b), which simpifies bookkeeping
at boundaries of compound rules.

W. Boland & C. Duris
Product Type Quadrature Formulas
BIT, 11, 139-158 (1971)

Created 2012-10-25 by Tom Loredo
"""

from numpy import sqrt, empty, array, concatenate

from quad import Quad, CompositeQuad


# Constants for Gauss-Legendre nodes:
rrt3 = 1/sqrt(3)  # reciprocal root 3
rt3_5 = sqrt(3./5)


class ProdQuad11(object):
    """
    Symmetric interpolatory (1,1) product quadrature rule.
    """
 
    def __init__(self, a, b, u_0, u_1, v_0=None, v_1=None,
                 f=None, fvals=None, ufunc=True):
        """
        Set up an interpolatory product quadrature rule over [a, b] using the
        specified nodes for the f(x) (u's) and g(x) (v's) factors.

        If the g(x) nodes are unspecified, use Gauss-Legendre nodes.

        If f or fvals is provided (f(x) at the u nodes), use it to build a
        rule for integrating g(x) specified by itself.

        When f is provided, ufunc indicates whether it is broadcastable.
        """
        self.a, self.b = a, b
        self.u_0, self.u_1 = u_0, u_1
        self.fnodes = array([u_0, u_1], float)
        if v_0 is None:
            if v_1 is None:  # use Gauss-Legendre nodes for g()
                mid = 0.5*(a + b)
                offset = 0.5*(b - a)*rrt3
                v_0 = mid - offset
                v_1 = mid + offset
            else:
                raise ValueError('Invalid v nodes!')
        else:
            self.v_0, self.v_1 = v_0, v_1
        self.gnodes = array([v_0, v_1], float)

        # Constants used in the weight matrix:
        self.ba = b - a
        self.ba2 = 0.5*(b**2 - a**2)
        self.ba3 = (b**3 - a**3) / 3.

        self.du = u_1 - u_0
        self.dv = v_1 - v_0

        # Calculate the weight matrix.
        self.fg_wts = empty((2,2))
        self.fg_wts[0,0] = self.wt_ij(u_0, u_1, v_0, v_1)
        self.fg_wts[0,1] = self.wt_ij(u_0, u_1, v_1, v_0)
        self.fg_wts[1,0] = self.wt_ij(u_1, u_0, v_0, v_1)
        self.fg_wts[1,1] = self.wt_ij(u_1, u_0, v_1, v_0)

        if f is None and fvals is None:
            self.fvals = None
        else:
            self.set_f(f, fvals, ufunc)

    def wt_ij(self, u, up, v, vp):
        """
        Calculate an element of the weight matrix.
        """
        return (self.ba3 - (up+vp)*self.ba2 + up*vp*self.ba) /\
               ((u - up)*(v - vp))

    def set_f(self, f=None, fvals=None, ufunc=True):
        """
        Calculate weights summed over tabulated values of f(x).
        """
        if fvals is None:
            if f is None:
                raise ValueError('Provide one of f or fvals!')
            if ufunc:
                self.fvals = f(self.fnodes)
            else:
                self.fvals = array( [f(self.fnodes[i]) for i in (1,2)] )
        else:
            if fvals is None:
                raise ValueError('Provide one of f or fvals!')
            self.fvals = fvals

        self.f = f
        self.gwts = empty(2)
        self.gwts[0] = self.fvals[0]*self.fg_wts[0,0] + self.fvals[1]*self.fg_wts[1,0]
        self.gwts[1] = self.fvals[0]*self.fg_wts[0,1] + self.fvals[1]*self.fg_wts[1,1]

    def quad_fg(self, f, g, ufunc=(True, True)):
        """
        Evaluate the quadrature rule using the functions f() and g().
        """
        if ufunc[0]:
            fvals = f(self.fnodes)
        else:
            fvals = array( [f(self.fnodes[i]) for i in (1,2)] )
        if ufunc[1]:
            gvals = g(self.gnodes)
        else:
            gvals = array( [g(self.gnodes[i]) for i in (1,2)] )
        return fvals[0]*self.fg_wts[0,0]*gvals[0] + fvals[0]*self.fg_wts[0,1]*gvals[1] + \
               fvals[1]*self.fg_wts[1,0]*gvals[0] + fvals[1]*self.fg_wts[1,1]*gvals[1]

    def quad_g(self, g, ufunc=True):
        """
        Evaluate the quadrature rule using the function g(), and previously
        specified f() values.
        """
        if self.fvals is None:
            raise ValueError('Unspecified f() values!')
        if ufunc:
            gvals = g(self.gnodes)
        else:
            gvals = array( [g(self.gnodes[i]) for i in (1,2)] )
        return self.gwts[0,0]*gvals[0] + self.gwts[1]*gvals[1]

    def quad_object(self, f=None, fvals=None, ufunc=True):
        """
        Return an object with the Quad interface interface expected by
        Composite quadtrature objects, using a specified f() (or a set of its
        values on the fnodes) to define a quadrature rule for g().
        """
        if f is not None or fvals is not None:  # otherwise assume prev. set
            self.set_f(f, fvals, ufunc)
        return Quad(self.a, self.b, self.gnodes, self.gwts)
