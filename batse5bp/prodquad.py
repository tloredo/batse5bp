""" 
Basic inner product quadrature algorithms, for integrating the product of a
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

The main reference guiding this implementation is:

W. Boland & C. Duris
Product Type Quadrature Formulas
BIT, 11, 139-158 (1971)

Created 2012-10-25 by Tom Loredo
"""

from numpy import sqrt, empty, ones_like, array, concatenate, arange
from numpy import sum, dot, prod, delete

from quad import Quad, CompositeQuad
from poly import roots2coefs, indef2def


# Constants for Gauss-Legendre nodes:
rrt3 = 1/sqrt(3)  # reciprocal root 3
rt3_5 = sqrt(3./5)


class ProdQuad1m(object):
    """
    Base class for interpolatory (1,m) product quadrature rules that calculate
    the integral of the product of two functions, f(x)*g(x).
    """
 
    def __init__(self, a, b, u_0, u_1, m):
        """
        Set up an interpolatory product quadrature rule over [a, b] using the
        specified nodes (u's) for f(x); subclasses specify g node behavior.

        m is the order for the g(x) dimension; there will be m+1 g nodes.
        """
        self.a, self.b = float(a), float(b)
        self.u_0, self.u_1 = u_0, u_1
        self.fnodes = array([u_0, u_1], float)
        self.nf = 2
        self.ng = m + 1

        self.fvals = None

    def wt_ij(self, u, up):
        """
        Calculate an element of the weight matrix.
        """
        raise NonImplementedError()

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
                self.fvals = array( [f(self.fnodes[i]) for i in range(self.nf)] )
        else:
            if fvals is None:
                raise ValueError('Provide one of f or fvals!')
            self.fvals = fvals

        self.f = f
        self.gwts = empty(self.ng)
        for i in range(self.ng):
            self.gwts[i] = sum(self.fvals*self.fg_wts[:,i])

    def quad_fg(self, f, g, ufunc=(True, True)):
        """
        Evaluate the quadrature rule using the functions f() and g().
        """
        if ufunc[0]:
            fvals = f(self.fnodes)
        else:
            fvals = array( [f(self.fnodes[i]) for i in range(self.nf)] )
        if ufunc[1]:
            gvals = g(self.gnodes)
        else:
            gvals = array( [g(self.gnodes[i]) for i in range(self.ng)] )
        return dot(fvals, dot(self.fg_wts, gvals))

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
            gvals = array( [g(self.gnodes[i]) for i in range(self.ng)] )
        return dot(self.gwts, gvals)

    def quad_object(self, f=None, fvals=None, ufunc=True):
        """
        Return an object with the Quad interface interface expected by
        Composite quadtrature objects, using a specified f() (or a set of its
        values on the fnodes) to define a quadrature rule for g().
        """
        if f is not None or fvals is not None:  # otherwise assume prev. set
            self.set_f(f, fvals, ufunc)
        return Quad(self.a, self.b, self.gnodes, self.gwts)


class ProdQuad11(ProdQuad1m):
    """
    Interpolatory (1,1) product quadrature rule.
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
        ProdQuad1m.__init__(self, a, b, u_0, u_1, 1)  # set up for m=1
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

        # Calculate the weight matrix.
        self.fg_wts = empty((2,2))
        self.fg_wts[0,0] = self.wt_ij(u_0, u_1, v_0, v_1)
        self.fg_wts[0,1] = self.wt_ij(u_0, u_1, v_1, v_0)
        self.fg_wts[1,0] = self.wt_ij(u_1, u_0, v_0, v_1)
        self.fg_wts[1,1] = self.wt_ij(u_1, u_0, v_1, v_0)

        if f is not None or fvals is not None:
            self.set_f(f, fvals, ufunc)

    def wt_ij(self, u, up, v, vp):
        """
        Calculate an element of the weight matrix.
        """
        return (self.ba3 - (up+vp)*self.ba2 + up*vp*self.ba) /\
               ((u - up)*(v - vp))


class ProdQuad12(ProdQuad1m):
    """
    Interpolatory (1,2) product quadrature rule.
    """
 
    def __init__(self, a, b, u_0, u_1, v_0=None, v_1=None, v_2=None,
                 f=None, fvals=None, ufunc=True):
        """
        Set up an interpolatory product quadrature rule over [a, b] using the
        specified nodes for the f(x) (u's) and g(x) (v's) factors.

        If the g(x) nodes are unspecified, use Gauss-Legendre nodes.

        If f or fvals is provided (f(x) at the u nodes), use it to build a
        rule for integrating g(x) specified by itself.

        When f is provided, ufunc indicates whether it is broadcastable.
        """
        ProdQuad1m.__init__(self, a, b, u_0, u_1, 2)  # set up for m=1
        if v_0 is None and v_1 is None and v_2 is None:  # use Gauss-Legendre nodes for g()
            mid = 0.5*(a + b)
            offset = 0.5*(b - a)*rt3_5
            v_0 = mid - offset
            v_1 = mid
            v_2 = mid + offset
        else:
            self.v_0, self.v_1, self.v_2 = v_0, v_1, v_2
        self.gnodes = array([v_0, v_1, v_2], float)

        # Constants used in the weight matrix:
        self.ba = b - a
        self.ba2 = 0.5*(b**2 - a**2)
        self.ba3 = (b**3 - a**3) / 3.
        self.ba4 = (b**4 - a**4) / 4.

        # Calculate the weight matrix.
        self.fg_wts = empty((2,3))
        self.fg_wts[0,0] = self.wt_ij(u_0, u_1, v_0, v_1, v_2)
        self.fg_wts[0,1] = self.wt_ij(u_0, u_1, v_1, v_0, v_2)
        self.fg_wts[0,2] = self.wt_ij(u_0, u_1, v_2, v_0, v_1)
        self.fg_wts[1,0] = self.wt_ij(u_1, u_0, v_0, v_1, v_2)
        self.fg_wts[1,1] = self.wt_ij(u_1, u_0, v_1, v_0, v_2)
        self.fg_wts[1,2] = self.wt_ij(u_1, u_0, v_2, v_0, v_1)

        if f is not None or fvals is not None:
            self.set_f(f, fvals, ufunc)

    def wt_ij(self, u, up, v, vp, vpp):
        """
        Calculate an element of the weight matrix.
        """
        return (self.ba4 - (up+vp+vpp)*self.ba3 \
                + (up*vp+up*vpp+vp*vpp)*self.ba2 - up*vp*vpp*self.ba) /\
               ((u - up)*(v - vp)*(v - vpp))


class ProdQuadRule(object):
    """
    An implementation of an interpolatory inner product quadrature rule for
    computing the integral of the product of two functions, f(x)*g(x),
    with support for alteration of the integration bounds for the rule without
    computing an entirely new rule from scratch.

    This implementation presumes the rule will be used to compute several
    quadratures (e.g., for different choices of f or g, or different intervals),
    and thus opts for simplicity over efficiency in initializing the rule.
    """

    # TODO:  Symmetric rules (with the same f and g nodes) can be computed
    # more quickly.  Support this, perhaps with a separate class?

    def __init__(self, f_nodes, g_nodes, a=None, b=None, f=None):
        """
        Define an interpolatory inner product quadrature rule for calculating
        the integral of the product of two functions, f(x)*g(x), over an
        interval [a,b].  The rule uses values of f() and g() at different
        sets of nodes.  Information is cached to enable defining new rules
        using the same f() nodes but different integration limits (and 
        different g() nodes).

        If limits are specified, a complete rule is constructed; in this
        case, methods are available for computing the quadrature.

        If `f` is provided, it is used as the first factor, f(x), defining an
        accelerated rule for quadrature when g() alone is subsequently 
        specified.  `f` may be either the (callable) function f(x), or a
        vector of values of f() evaluated over f_nodes.  It is used only if
        [a,b] is specified.  f() may also be set later via the set_f() method.
        """
        # The # of nodes minus 1 is the degree of the Lagrange polynomial used
        # for interpolation.  Typically this is the degree (separate for f & g)
        # for which quadrature of polynomials is exact.  If the nodes are
        # Gaussian quadrature nodes, the result will be exact for polynomials of
        # higher degree (with appropriate weight function).
        self.n_f = len(f_nodes)
        self.f_deg = self.n_f - 1
        self.n_g = len(g_nodes)
        self.g_deg = self.n_g - 1
        self.f_nodes = array(f_nodes, dtype=float)
        self.g_nodes = array(g_nodes, dtype=float)
        if a is None or b is None:
            if a is not None or b is not None:
                raise ValueError('Must specify *both* a and b!')
            self.a, self.b = None, None
            bounded = False
            self.fvals = None
        else:
            self.a, self.b = a, b
            bounded = True
            if f:
                if callable(f):
                    self.fvals = f(f_nodes)
                else:
                    self.fvals = f.copy()
            else:
                self.fvals = None

        # The rule is specified by a (n_f+1)x(n_g+1) matrix of weights.
        # The denominators do not depend on [a,b] and so are cached for
        # use when [a,b] is changed.
        self.fg_wts = empty((self.n_f,self.n_g), dtype=float)
        self.f_denom = ones_like(self.f_nodes)
        for i in xrange(self.n_f):
            self.f_denom[i] = prod(self.f_nodes[i] - delete(self.f_nodes, i))
            # fd = 1
            # for k in xrange(self.n_f):
            #     if k != i:
            #         fd *= (self.f_nodes[i] - self.f_nodes[k])
        self.g_denom = ones_like(self.g_nodes)
        for j in xrange(self.n_g):
            self.g_denom[j] = prod(self.g_nodes[j] - delete(self.g_nodes, j))

        # The weights are integrals of products of the Lagrange basis
        # polynomials over [a,b].  Store arrays defining the indefinite
        # integrals (via polynomial coefficients) so weights for new [a,b] can
        # be quickly calculated.
        # Each indefinite integral is a f_deg*g_deg+1 degree polynomial without
        # a constant term, thus with f_deg+g_deg+1 coefficients.
        n_c = self.f_deg+self.g_deg+1
        self.indefs = empty((self.n_f, self.n_g, n_c), dtype=float)

        # powers will hold divisors converting polynomial coefs to indefinite
        # integral coefs, i.e., the powers of the integrated monomial terms.
        powers = arange(1, n_c+1, dtype=float)
        f_roots = empty(self.f_deg)
        g_roots = empty(self.g_deg)
        roots = empty(self.f_deg+self.g_deg)
        j0 = self.n_f - 1  # index for the first g_roots entry in roots
        # print
        # print self.n_f, self.n_g
        for i in xrange(self.n_f):
            if i>0:
                roots[0:i] = self.f_nodes[0:i]
            if i<self.n_f-1:
                # Start at i since we're skipping entry i.
                roots[i:j0] = self.f_nodes[i+1:]
            for j in xrange(self.n_g):
                if j>0:
                    roots[j0:j0+j] = self.g_nodes[0:j]
                if j<self.n_g-1:
                    # Start at j0+j since we're skipping entry j.
                    roots[j0+j:] = self.g_nodes[j+1:]
                # print i, j, roots
                coefs = roots2coefs(roots)
                self.indefs[i,j,:] = coefs/powers

        if bounded:
            self.set_wts(a, b)

        # *** This duplicates some initialization done above.
        if self.fvals is not None:
            self.set_f(fvals=self.fvals)

    def set_wts(self, a, b):
        """
        Using the integration range [a,b], calculate the quadrature weight
        matrix.
        """
        for i in xrange(self.n_f):
            for j in xrange(self.n_g):
                # TODO:  This may be wasteful, repeatedly calculating powers of
                # a and b.
                self.fg_wts[i,j] = indef2def(self.indefs[i,j,:], a, b) / \
                                (self.f_denom[i]*self.g_denom[j])

    def get_fg_wts(self, a, b):
        """
        Using the integration range [a,b], calculate the quadrature weight
        matrix and return it.  This does not change the value of any
        previously-set weight matrix stored in the instance, and thus has
        no effect on future calls of the quad_fg() and quad_g() methods.
        """
        fg_wts = empty((self.n_f,self.n_g), dtype=float)
        for i in xrange(self.n_f):
            for j in xrange(self.n_g):
                # TODO:  This may be wasteful, repeatedly calculating powers of
                # a and b.
                fg_wts[i,j] = indef2def(self.indefs[i,j,:], a, b) / \
                              (self.f_denom[i]*self.g_denom[j])
        return fg_wts

    def set_f(self, f=None, fvals=None, ufunc=True):
        """
        Setup a weight matrix using tabulated values of f(x), for subsequent
        use for quadratures with specified g(x) functions.
        """
        if fvals is None:
            if f is None:
                raise ValueError('Provide one of f or fvals!')
            if ufunc:
                self.fvals = f(self.f_nodes)
            else:
                self.fvals = array( [f(self.f_nodes[i]) for i in range(self.n_f)] )
        else:
            if fvals is None:
                raise ValueError('Provide one of f or fvals!')
            self.fvals = fvals

        self.f = f
        self.g_wts = empty(self.n_g)
        for j in range(self.n_g):
            self.g_wts[j] = sum(self.fvals*self.fg_wts[:,j])

    def quad_fg(self, f, g, ufunc=(True, True)):
        """
        Evaluate the quadrature rule using the functions f() and g().
        """
        if ufunc[0]:
            fvals = f(self.f_nodes)
        else:
            fvals = array( [f(self.f_nodes[i]) for i in range(self.n_f)] )
        if ufunc[1]:
            gvals = g(self.g_nodes)
        else:
            gvals = array( [g(self.g_nodes[j]) for j in range(self.n_g)] )
        return dot(fvals, dot(self.fg_wts, gvals))

    def quad_g(self, g, ufunc=True):
        """
        Evaluate the quadrature rule using the function g(), and previously
        specified f() values.
        """
        if self.fvals is None:
            raise ValueError('Unspecified f() values!')
        if ufunc:
            gvals = g(self.g_nodes)
        else:
            gvals = array( [g(self.g_nodes[j]) for j in range(self.n_g)] )
        return dot(self.gwts, gvals)

    def quad_object(self, f=None, fvals=None, ufunc=True):
        """
        Return an object with the Quad interface interface expected by
        Composite quadtrature objects, using a specified f() (or a set of its
        values on the fnodes) to define a quadrature rule for g().
        """
        if f is not None or fvals is not None:  # otherwise assume prev. set
            self.set_f(f, fvals, ufunc)
        return Quad(self.a, self.b, self.g_nodes, self.g_wts)
