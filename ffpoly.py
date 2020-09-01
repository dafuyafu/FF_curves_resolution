from sympy.polys.polytools import Poly, poly, LM, LT
from sympy.core.numbers import Integer
from sympy.core.expr import Expr

class SFF:
    """
        represents a splitted finite field of some polynomials over a finite field of which modulus is 'mod'.
    """

    def __init__(self, rel, mod):
        """

        Instance variables:
            * rel_list: list of relational equations
                        of which element is a dict {'var': variable, 'rep': equation}
            * mod: characteristic number

        Example:
            a_1 = symbols('a_1') 
            rel = {'var': a_1, 'rep': a_1 ** 2 - 3, 'deg': 2}
            rel_list = [rel, rel2, ...]

        """
        if isinstance(rel, Expr):
            self.rel_list = [{'var': poly(rel).gens[0], 'rep': rel, 'deg': poly(rel).degree()}]
        elif isinstance(rel, list):
            self.rel_list = []
            for _p in rel:
                self.rel_list.append({'var': poly(_p).gens[0], 'rep': _p, 'deg': poly(rel).degree()})
        else: 
            raise ValueError("first argument needs to be an Expr instance or list.")
        self.mod = mod
        _var_list = []
        for _rel in self.rel_list:
            if not _rel['var'] in _var_list:
                _var_list.append(_rel['var'])
        self.gens = _var_list

    def as_FF(self):
        """ Not implemented yet """
        """ output FF(n^r) """
        raise NotImplementedError("unchi!")

    def as_sympy_FF(self):
        return "FF(" + str(self.mod) + ")"

    def as_SFF(self):
        return "SFF(" + str(self.mod) + ") with " + str(self.rel_list)

    def ext_deg(self):
        """ Not implemented yet """
        raise NotImplementedError("unchi!")

    def rel_deg(self, **args):
        if not args:
            return poly(self.rel_list[0]['rep']).degree()
        elif 'var' in args:
            for rel in self.rel_list:
                if args['var'] == rel['var']:
                    return poly(rel['rep']).degree()
            raise ValueError("doesn't have the variable")
        elif 'index' in args:
            return poly(self.rel_list[args['index']]['rep']).degree()
        else:
            raise ValueError("rel_deg() doesn't the argument option")

    def rel_append(self, rel):
        """ Validations are not implemented. """
        if isinstance(rel, Expr):
            self.rel_list.append({'var': poly(rel).gens[0], 'rep': rel, 'deg': poly(rel).degree()})
        elif isinstance(rel, list):
            for _rel in rel:
                self.rel_list.append({'var': poly(_rel).gens[0], 'rep': _rel, 'deg': poly(_rel).degree()})
        else:
            pass

    def point_list(self):
        raise NotImplementedError()

    def point(self, i):
        


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
        if isinstance(rep, Poly):
            self.rep = rep.as_expr()
        else:
            self.rep = rep
        self.dom = dom
        self.var = poly(rep).gens

        if isinstance(rep, Integer):
            self.is_int = True
        else:
            self.is_int = False

    def __str__(self):
        return "FFPoly(" + str(self.rep) + ", modulus=" + str(self.dom.mod) + ", rel: " + str(self.dom.rel_list) + ")"

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

    @classmethod
    def ffpoly(rep, dom):
        return FFPoly(rep, dom)

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
            return ffpoly(_add.as_expr(), f.dom)
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
            return ffpoly(_sub.as_expr(), f.dom)
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
            return ffpoly(_mul.as_expr(), f.dom)._reduce()
        else:
            raise ValueError("cannot add polynomials over different FFs")

    def _reduce(self):
        _ffpoly = self
        for rel in self.dom.rel:
            _ffpoly = _ffpoly._simple_reduce(rel)
        return _ffpoly

    def simple_reduce(self, rel):
        return self._simple_reduce(rel)

    def _simple_reduce(self, rel):
        """
        Divide f by its relation polynomial and calculate the residue
        when f.rel is monovariate i.e. considering over simple extention field.

        Examples
        ========
        >>> f = SFFPoly(a ** 2 * x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> f._simple_reduce()
        SFFPoly(-2 * x + a * y, a ** 2 - 2, modulus=5)
        """
        _lm = LM(poly(rel['rep']))
        _sub = _lm - rel['rep']

        if self.degree(rel['var']) < poly(rel['rep']).degree(rel['var']):
            """ _rel doesn't devide _p """
            return self
        else:
            return self._simple_reduce_rec(_lm, _sub, rel['var'])

    def _simple_reduce_rec(self, lm, sub, var):
        if self.degree(var) == poly(lm).degree(var):
            return ffpoly(self.poly.subs({lm: sub}), dom=self.dom)
        else:
            _p = self._simple_reduce_rec(lm * var, sub * var, var).rep.subs({lm: sub})
            return ffpoly(_p, dom=self.dom)

    def solve(self):


    def solveabs(self):
        """ solve self over its algebraic closure """
        raise NotImplementedError

class P2Point:
    """
    represents a point of P2 over finite field.
    """

    def __init__(self, a, b, c, dom):
        if a == 0 and b == 0 and c == 0:
            raise ValueError("(0,0,0) is not a P2 point")
        x,y,z = symbols('x,y,z')
        self.cod = {x: a, y: b, z: c}
        self.dom = dom

    def __str__(self):
        return str(self.coordinate) + ' (rel: ' + str(self.rel) + ')'

    def __eq__(p, q):
        if 



