"""
Generic basic and composite quadrature rule objects.

Created 2012-10-28 by Tom Loredo; stolen from the Inference package
"""

from numpy import concatenate


class Quad:
    """
    A simple quadrature rule
    """
    
    def __init__(self, *args):
        """
        Define a quadrature rule from its nodes and weights.
        
        For closed rules (where the nodes include the boundaries), the
        signature may be:
        
            QuadRule(nodes, wts)
        
        where `nodes` and `wts` are arrays containing the (ordered)
        nodes and weights for the rule.
        
        For open rules (and optionally for closed rules), the signature is
        
            QuadRule(l, u, nodes, wts)

        where `a` and `b` are the integration limits, and `nodes` and `wts` 
        are arrays containing the (ordered) nodes and weights for the rule.
        """
        if len(args) == 2:
            self.type = 'closed'
            self.nodes = args[0]
            self.wts = args[1]
            self.l = self.nodes[0]
            self.u = self.nodes[-1]
        elif len(args) == 4:
            self.l = args[0]
            self.u = args[1]
            self.nodes = args[2]
            self.wts = args[3]
            if self.l==self.nodes[0] and self.u==self.nodes[-1]:
                self.type = 'closed'
            else:
                self.type = 'open'
        else:
            raise ValueError('Specify (nodes, wts) or (a, b, nodes, wts)!')
        self.npts = len(self.nodes)
        if len(self.wts) != self.npts:
            raise ValueError('Mismatch in length of nodes and wts!')

    def quad(self, f):
        """
        Evaluate the quadrature given a function or array of function values.
        
        If f is callable, find the array of values f(self.nodes) and
        evaluate the quadrature.  Otherwise, f must be a vector of function
        values at the nodes, used to evaluate the quadrature.
        """
        if callable(f):
            fvals = f(self.nodes)
        elif len(f) == self.npts:
            fvals = f
        else:
            raise ValueError('Argument must be callable or array of values!')
        return sum(self.wts*fvals)


class CompositeQuad:
    """
    Composite quadrature rule built over contiguous intervals
    """

    @staticmethod
    def isquad(obj):
        """
        Return True if obj has a Quad interface.
        """
        if isinstance(obj, Quad):
            return True
        try:
            assert hasattr(obj, 'l')
            assert hasattr(obj, 'u')
            assert hasattr(obj, 'npts')
            assert hasattr(obj, 'nodes')
            assert hasattr(obj, 'wts')
            return True
        except AssertionError:
            return False

    def __init__(self, *args):
        """
        Define a composite quadrature rule from a collection of rules.
        
        Each argument should describe a basic quadrature rule.  There are
        three possible formats:
        
          * The argument may be a QuadRule instance.
        
          * For a closed rule, the argument may be a 2-tuple of the
            form (nodes, wts), where `nodes` and `wts` are arrays of
            values for the nodes and weights of the rule.
          
          * For an open or closed rule, the argument may be a 4-tuple of the
            form (a, b, nodes, wts), where `a` and `b` are the limits of
            integration, and `nodes` and `wts` are as above.
        
        The rules must be contiguous.
        """
        # Collect all the rules as QuadRule instances.
        self.rules = []
        last = None
        for i, arg in enumerate(args):
            if CompositeQuad.isquad(arg):
                if last:
                    if last.u != arg.l:
                        raise ValueError(\
                            'Rule %i not contiguous with previous!' % (i+1))
                self.rules.append(arg)
                last = arg
            else:
                rule = Quad(*arg)
                self.rules.append(rule)
                last = rule
        # Make an array of all the distinct nodes.
        prev = self.rules[0]
        distinct = [ prev.nodes ]
        self.starts = [0]  # keeps track of starting indices for rules
        self.npts = prev.npts
        for rule in self.rules[1:]:
            if prev.nodes[-1] != rule.nodes[0]:
                distinct.append(rule.nodes)
                self.starts.append(self.npts)
                self.npts += rule.npts
            else:
                distinct.append(rule.nodes[1:])
                self.starts.append(self.npts-1)
                self.npts += rule.npts - 1
            prev = rule
        self.nodes = concatenate(distinct)
        self.factor = None
        self.factors = None
        self.l = self.rules[0].l
        self.u = self.rules[-1].u

    def quad(self, f):
        """
        Evaluate the quadrature of the callable f.
        """
        fvals = f(self.nodes)
        self.fvals = fvals
        if self.factor:
            self.ivals = self.factor * self.fvals
        else:
            self.ivals = self.fvals
        result = 0.
        # Go through the rules, passing the function evaluations.
        for rule, start in zip(self.rules, self.starts):
            result += rule.quad(self.ivals[start:start+rule.npts])
        return result

