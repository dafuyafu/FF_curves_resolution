from sympy.polys.polytools import Poly, LM, LT

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
        if len(rel.gens) != 1:
            raise ValueError("relational polynomial need to be monovariate")
        if LT(rel) != LM(rel):
            raise ValueError("relational polynomial need to be monic")
        return super().__new__(cls, rep, *gens, **args)

    def __init__(self, rep, rel, *gens, **args):
        if rel == None:
            raise ValueError("need a relation polynomial")
        else:
            self.rel = rel
            self.rel_var = rel.gens[0]
            self.poly = Poly(rep, *gens, **args)

    def __str__(self):
        print("SFF" + self.poly.__str__() + " on " + str(self.rel))

    def add(f,g):
        """
        Add two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = SFF(x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = SFF(x + 4 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f + g
        SFF(2 * x, x, modulus=5)
        """
        if f.domain != g.domain:
            raise ValueError("cannot add polynomials over different FFs")
        if f.rel == g.rel:
            _add = f.poly + g.poly
            return SFFPoly(_add, f.rel, domain=f.domain)._simple_reduce()
        else:
            raise ValueError("cannot add polynomials over different SFFs")

    def sub(f,g):
        """
        Add two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = SFF(x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = SFF(x + 4 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f + g
        SFF(2 * x, x, modulus=5)
        """
        if f.domain != g.domain:
            raise ValueError("cannot add polynomials over different FFs")
        if f.rel == g.rel:
            _sub = f.poly - g.poly
            return SFFPoly(_sub, f.rel, domain=f.domain)._simple_reduce()
        else:
            raise ValueError("cannot add polynomials over different SFFs")

    def mul(f,g):
        """
        Add two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = SFF(x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = SFF(x + 4 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f + g
        SFF(2 * x, x, modulus=5)
        """
        if f.domain != g.domain:
            raise ValueError("cannot add polynomials over different FFs")
        if f.rel == g.rel:
            _mul = f.poly * g.poly
            return SFFPoly(_mul, f.rel, domain=f.domain)._simple_reduce()
        else:
            raise ValueError("cannot add polynomials over different SFFs")

    def _simple_reduce(self):
        """
        Divide f by its relation polynomial and calculate the residue
        when f.rel is monovariate i.e. considering over simple extention field.

        Examples
        ========
        >>> f = SFF(a ** 2 * x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> f._simple_reduce()
        SFF(-2 * x + a * y, a ** 2 - 2, modulus=5)
        """
        if self.rel_var not in self.poly.gens:
            return self

        _lm = LM(self.rel)
        _sub = _lm - self.rel.as_expr()

        if self.degree(self.rel_var) < self.rel.degree(self.rel_var):
            """ _rel doesn't devide _p """
            return self
        else:
            return self._simple_reduce_rec(_lm, _sub)

    def _simple_reduce_rec(self, lm, sub):
        if self.degree(self.rel_var) == Poly(lm).degree(self.rel_var):
            return SFFPoly(self.poly.subs({lm: sub}), self.rel, domain=self.domain)
        else:
            _p = self._simple_reduce_rec(lm * self.rel_var, sub * self.rel_var).poly.subs({lm: sub})
            return SFFPoly(_p, self.rel, domain=self.domain)
