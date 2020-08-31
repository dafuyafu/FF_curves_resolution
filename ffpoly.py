from sympy.polys.polytools import poly, LM, LT
from sympy.core.numbers import Integer

class SFF:
    def __init__(self, rel, mod):
        if isinstance(rel, Expr):
            self.rel = [rel]
        elif isinstance(rel, list):
            self.rel = rel
        else: 
            raise ValueError("first argument needs to be an Expr instance or list.")
        self.mod = mod

    def as_FF(self):
        """Not implemented yet"""
        """output FF(n^r)"""

    def as_sympy_FF(self):
        return "FF(" + str(self.mod) + ")"

    def as_SFF(self):
        return "SFF(" + str(self.mod) + ") with " + str(self.rel)

    def rel_var(self)
        _list = []
        for _p in self.rel:
            for _v in poly(_p).gens:
                if not _v in _list:
                    _list.append(_v)
        return tuple(_list)

    def degree_rel(self, i):
        return poly(self.rel[i]).degree()

class FFPoly:
    """ 
    represents an element of polynomial ring

    Examples
    ========

    >>> from sympy.polys.polytools import Poly
    >>> from sympy.cores.symbol import symbol

    Create a FFPoly instance:

    >>> x,y,z,a = symbols('x,y,z,a')
    >>> _f = y ** 2 * z ** 3 + 2 * x ** 5 + 2 * x ** 3 * z ** 2 - 2 * x * z ** 4
    >>> f = FFPoly(_f, a ** 2 - 2, domain='FF(5)')

    """

    def __init__(self, rep, dom):
        """
            ## arguments:

            * rep: expression of the polynomial
            * dom: domain field which is a SFF instance

        """
        self.rep = rep
        self.dom = dom
        self.var = poly(f).gens

        if isinstace(rep, Integer)
            self.is_int = True
        else:
            self.is_int = False

    def __str__(self):
        print("FFPoly(" + str(self.rep) + ", modulus=" + str(self.dom.mod) + ", rel= " + str(self.dom.rel) + ")")

    def as_expr(self):
        return self.rep

    def as_poly(self):
        if self.is_int:
            raise ValueError("this is an integer.")
        else:
            return poly(self.rep, domain=self.as_sympy_FF)

    def degree(self, *gens):
        if self.rep == 0:
            return "-oo"
        elif self.is_int:
            return 0
        elif len(gens) == 0:
            return poly(self.rep).degree()
        elif len(gens) == 1:
            return poly(self.rep).degree(gens[0])
        else:
            raise ValueError("need only one or zero argument")

    def __add__(f,g):
        """
        Add two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = FFPoly(x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = FFPoly(x + 4 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f + g
        SFFPoly(2 * x, x, modulus=5)
        """
        if f.dom == g.dom:
            _add = f.as_poly() + g.as_poly()
            return FFPoly(_add, f.dom)
        else:
            raise ValueError("cannot add polynomials over different FFs")

    def sub(f,g):
        """
        subtract two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = SFFPoly(2 * x - a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = SFFPoly(x + 4 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f - g
        SFFPoly(x, x, a, modulus=5)
        """
        if f.dom == g.dom:
            _sub = f.as_poly() - g.as_poly()
            return FFPoly(_sub, f.dom)
        else:
            raise ValueError("cannot add polynomials over different FFs")

    def mul(f,g):
        """
        multiple two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = SFFPoly(x - a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = SFFPoly(x - 2 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f * g
        SFFPoly(x ** 2 + 2 * a * x * y - y ** 2, x, modulus=5)
        """
        if f.dom == g.dom:
            _mul = f.poly * g.poly
            return FFPoly(_mul, f.dom)._reduce()
        else:
            raise ValueError("cannot add polynomials over different SFFs")

    def _reduce(self):

    def simple_reduce(self, var):
        return self._simple_reduce(var)

    def _simple_reduce(self, var):
        """
        Divide f by its relation polynomial and calculate the residue
        when f.rel is monovariate i.e. considering over simple extention field.

        Examples
        ========
        >>> f = SFFPoly(a ** 2 * x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> f._simple_reduce()
        SFFPoly(-2 * x + a * y, a ** 2 - 2, modulus=5)
        """

        _lm = LM(poly(self.rel))
        _sub = _lm - self.rel

        if self.degree(var) < self.rel.degree(var):
            """ _rel doesn't devide _p """
            return self
        else:
            return self._simple_reduce_rec(_lm, _sub)

    def _simple_reduce_rec(self, lm, sub):
        if self.degree(self.rel_var) == Poly(lm).degree(self.rel_var):
            return FFPoly(self.poly.subs({lm: sub}), self.rel, domain=self.domain)
        else:
            _p = self._simple_reduce_rec(lm * self.rel_var, sub * self.rel_var).poly.subs({lm: sub})
            return FFPoly(_p, self.rel, domain=self.domain)

class P2Point:
    """
    represents a point of P2 over finite field.
    """

    def __init__(self, a, b, c, rel):
        if a == 0 and b == 0 and c == 0:
            raise ValueError("(0,0,0) is not a PP2 point")
        x,y,z = symbols('x,y,z')
        self.coordinate = {x: a, y: b, z: c}
        self.rel = rel

    def __str__(self):
        print(str(self.coordinate) + ' (rel: ' + str(self.rel) + ')')