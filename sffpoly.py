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

    # def __str__(self):
    #     print(self.poly.__str__() + 'on SFF(%s, %s)', (self.poly.modulus, self.poly.rel))

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
        if f.rel == g.rel:
            _sum = f.poly + g.poly
            return SFFPoly(_sum, f.rel)
        else:
            raise ValueError("cannot add polynomials over different SFF")

    def _reduce(self):
        """
        Divide f by its relation polynomial and calculate the residue.

        Examples
        ========
        >>> f = SFF(a ** 2 * x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> f._reduce()
        SFF(-2 * x + a * y, a ** 2 - 2, modulus=5)
        """
        if self.rel_var not in self.poly.gens:
            return self

        _p = self.poly
        _lm = Poly(LM(self.rel), domain=self.domain)
        _rel = _lm - self.rel
        while(_p.degree(self.rel_var) >= _lm.degree()):
            _tlm = _lm
            _trel = _rel
            while(_p.degree(self.rel_var) > _tlm.degree()):
                _tlm *= self.rel_var
                _trel *= self.rel_var
            _p = Poly(_poly.subs({_tlm: _rel}), domain=self.domain)

        return SFFPoly(_poly, self.rel)

