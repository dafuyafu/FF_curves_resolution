from sympy.polys.polytools import Poly

class SFFPoly(Poly):
    """ 
    represents a polynomial over a splitting finite field.

    Examples
    ========

    >>> from sympy.polys.polytools import Poly
    >>> from sympy.cores.symbol import symbol

    Create a SFFPoly instance:

    >>> x,y,z,a = symbols('x,y,z,a')
    >>> _f = y ** 2 * z ** 3 + 2 * x ** 5 + 2 * x ** 3 * z ** 2 - 2 * x * z ** 4
    >>> f = SFFPoly(_f, a ** 2 - 2, domain='FF(5)')

    """

    def __new__(cls, rep, rel, *gens, **args):
        """Create a new polynomial instance via Poly.__new__()"""
        return super().__new__(cls, rep, *gens, **args)

    def __init__(self, rep, rel, *gens, **args):
        if rel == None:
            raise ValueError("Need a relational polynomial.")
        else:
            self.rel = rel

    def is_the_same_domain_as(f,g):
        return f.rel == g.rel

    def add(f,g):
        """
        Add two polynomials ``f`` and ``g`` and divide it with rel polynomial

        Examples
        ========

        >>> f = SFF(x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = SFF(x + 4 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f + g
        SFF(2 * x, domain='FF(5)')
        """



