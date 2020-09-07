from itertools import product

from sympy.core.numbers import ilcm, Integer
from sympy.core.expr import Expr
from sympy.core.function import diff, expand
from sympy.ntheory.primetest import isprime
from sympy.polys.polytools import factor_list, LC, Poly, poly, resultant

class SFF:
    """
        represents a splitted finite field of some polynomials
        over a finite field of which modulus is 'mod'.
    """

    def __init__(self, rel, mod):
        """
        Instance variables:
            * rel_list: list of relational equations
                        of which element is a dict {'var': variable, 'rep': equation}
            * mod: characteristic number
            * var_list: list of variables of relational equations
            * exdeg: extension degree
            * num: number of elements
            * gens: list of generators as vector space

        Example:
            a_1 = symbols('a_1') 
            rel = {'var': a_1, 'rep': a_1 ** 2 - 3, 'deg': 2, 'uni': True}

            self.rel_list = [rel_1, rel_2, ..., rel_n]
            self.mod = p
            self.var_list = [a_1, a_2, ..., a_m]
            self.exdeg = e
            self.num = p ** e
            self.gens = [g_1, g_2, ..., g_(p ** e)]
        """
        self.rel_list = []
        _dom = 'FF(' + str(mod) + ')'
        if rel == 0:
            self.is_prime = True
        elif isinstance(rel, Expr):
            self.is_prime = False
            if not LC(rel.as_poly()) == 1:
                rel = poly(rel * pow(LC(rel.as_poly()), mod - 2, mod), domain=_dom).as_expr()
            if len(poly(rel).gens) == 1:
                self.rel_list.append({'var': poly(rel).gens[0],
                                      'rep': rel.as_expr(), 
                                      'deg': poly(rel).degree(), 
                                      'is_uni': True})
            else:
                self.rel_list.append({'var': poly(rel).gens[0], 
                                      'rep': rel.as_expr(), 
                                      'deg': poly(rel).degree(), 
                                      'is_uni': False})
        elif isinstance(rel, list):
            self.is_prime = False
            if rel == []:
                self.rel_list.append({'var': 1, 
                                      'rep': 0, 
                                      'deg': "-oo", 
                                      'is_uni': True})
            for _p in rel:
                if not LC(_p.as_poly()) == 1:
                    _p = poly(_p * pow(LC(_p.as_poly()), mod - 2, mod), domain=_dom).as_expr()
                if len(poly(_p).gens) == 1:
                    self.rel_list.append({'var': poly(_p).gens[0], 
                                          'rep': _p.as_expr(), 
                                          'deg': poly(_p).degree(), 
                                          'is_uni': True})
                else:
                    self.rel_list.append({'var': poly(_p).gens[0], 
                                          'rep': _p.as_expr(), 
                                          'deg': poly(_p).degree(), 
                                          'is_uni': False})
        else: 
            raise ValueError("first argument needs to be 0, an Expr instance or list.")

        if isprime(mod):
            self.mod = mod
        else:
            raise ValueError("modulus needs to be a prime number")

        if self.is_prime:
            self.var_list = []
            self.exdeg = 1
            self.num = self.mod
            self.gens = [1]
        else:
            _var_list, _degs, _gens = [], [], []
            for _rel in self.rel_list:
                _degs.append(_rel['deg'])
                if not _rel['var'] in _var_list and _rel['is_uni']:
                    _var_list.append(_rel['var'])
            self.var_list = _var_list
            self.exdeg = ilcm(1, *_degs)
            self.num = mod ** self.exdeg
            """
                Note:
                Below codes include errors.
                When degree of two extension variables are each 2 and 4 then returns 8 generators,
                however acutually there exist only 4 generators.
            """
            _expr = 1
            for _var in self.var_list:
                _deg = next(item for item in self.rel_list if item['var'] == _var and item['is_uni'])['deg']
                _expr_1 = 0
                for i in range(_deg):
                    _expr_1 += _var ** i
                _expr *= _expr_1
            _tuple = expand(_expr).as_coeff_add()
            self.gens = tuple([_tuple[0]] + [item for item in _tuple[1]])

    def __str__(self):
        return self.as_SFF()

    def __repr__(self):
        return self.as_SFF()

    @classmethod
    def sffpoly(cls, rel, mod):
        return sffpoly(cls, rel, mod)

    def as_FF(self):
        return "FF(%s ** %s)" % (self.mod, self.exdeg)

    def as_sympy_FF(self):
        return "FF(%s)" % self.mod

    def as_SFF(self):
        if self.is_prime:
            return self.as_sympy_FF()
        else:
            return "SFF(%s ** %s) splitting %s" % (self.mod, self.exdeg, [rel['rep'] for rel in self.rel_list if rel['is_uni']])

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
            raise ValueError("rel_deg() doesn't have the argument option")

    def rel_append(self, rel):
        """ Validations are not implemented. """
        """ 
            CAUTION!!
            DO NOT USE THIS METHOD
        """
        if isinstance(rel, Expr):
            self.rel_list.append({'var': poly(rel).gens[0], 'rep': rel, 'deg': poly(rel).degree()})
        elif isinstance(rel, list):
            for _rel in rel:
                self.rel_list.append({'var': poly(_rel).gens[0], 'rep': _rel, 'deg': poly(_rel).degree()})
        else:
            pass

    def element(self, i):
        if i > self.num:
            raise ValueError
        _rep = 0
        for d in range(self.exdeg)[::-1]:
            _rep += int(i / self.mod ** d) * self.gens[d]
            i %= self.mod ** d
        return _rep

    def element_itr(self):
        for i in range(self.num):
            yield self.element(i)

    def elements_list(self):
        return [self.element(i) for i in range(self.num)]

class SFFPoly:
    """ 
    represents an element of polynomial ring

    Examples
    ========

    >>> from sympy.polys.polytools import Poly
    >>> from sympy.cores.symbol import symbol

    Create a sffpoly instance:

    >>> x,y,z,a = symbols('x,y,z,a')
    >>> _f = y ** 2 * z ** 3 + 2 * x ** 5 + 2 * x ** 3 * z ** 2 - 2 * x * z ** 4
    >>> f = sffpoly(_f, a ** 2 - 2, domain='FF(5)')

    """

    def __init__(self, rep, dom):
        """
            Instance variables:
            * rep: expression of the polynomial
            * dom: domain field which is a SFF instance
            * var: list of variables of rep
                   if rep is an integer then var == []
            * is_uni: boolean value indicating if rep is univariate or not 
            * is_int: boolean value indicating if rep is an integer or not

        """
        
        if isinstance(rep, Integer):
            self.rep = rep
            self.var = []
            self.is_int = True
        else:
            self.rep = reduce(rep.as_expr(), dom)
            self.var = [v for v in poly(self.rep).gens if not v in dom.var_list] 
            self.is_int = False
        self.dom = dom
        if len(self.var) == 1:
            self.is_uni = True
        else:
            self.is_uni = False

    def __repr__(self):
        return "SFFPoly(%s, %s)" % (self.rep, self.dom.as_SFF())

    def __add__(f,g):
        """
        Add two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = sffpoly(x + a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = sffpoly(x + 4 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f + g
        Ssffpoly(2 * x, x, modulus=5)
        """
        if f.dom == g.dom:
            _add = poly(reduce(expand(f.rep + g.rep), f.dom), domain=f.dom.as_sympy_FF()).as_expr()
            return sffpoly(_add, f.dom)
        else:
            raise ValueError("argument sffpolys have different domains")

    def __sub__(f,g):
        """
        subtract two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = Ssffpoly(2 * x - a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = Ssffpoly(x + 4 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f - g
        Ssffpoly(x, x, a, modulus=5)
        """
        if f.dom == g.dom:
            _sub = poly(reduce(expand(f.rep - g.rep), f.dom), domain=f.dom.as_sympy_FF()).as_expr()
            return sffpoly(_sub, f.dom)
        else:
            raise ValueError("argument sffpolys have different domains")

    def __mul__(f,g):
        """
        multiple two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = Ssffpoly(x - a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = Ssffpoly(x - 2 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f * g
        Ssffpoly(x ** 2 + 2 * a * x * y - y ** 2, x, modulus=5)
        """
        if f.dom == g.dom:
            _mul = poly(reduce(expand(f.rep * g.rep), f.dom), domain=f.dom.as_sympy_FF()).as_expr()
            return sffpoly(_mul, f.dom)
        else:
            raise ValueError("argument sffpolys have different domains")

    def __eq__(f,g):
        if f.dom == g.dom and f.rep == g.rep:
            return True
        else:
            return False

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

    def subs(self, point):
        if self.is_int:
            raise TypeError("this is a modular integer")
        if point == {}:
            raise ValueError
        _rep = self.rep
        for k, v in point.items():
            _rep = _rep.subs({k: v})
        return _rep

    def subs_as_sffpoly(self, **args):
        return sffpoly(self.subs(args), self.dom)
        
    def solve_abs(self):
        """ solve self over its algebraic closure """
        raise NotImplementedError

    def simple_solve(self):
        _sol = []
        for point in product(self.dom.element_itr(), repeat=len(self.var)):
            point_dict = dict((_var, point[self.var.index(_var)]) for _var in self.var)
            if reduce(self.subs(point_dict), self.dom) == 0:
                _sol.append(point_dict)
        return _sol

    def diff(self, var):
        return reduce(diff(self.rep, var, self.dom))

    def sing(self):
        """ find singular locus of self """
        pass

class P2Point:
    """
    represents a point of P2 over finite field.
    """

    def __init__(self, cod, dom):
        x,y,z = symbols('x y z')
        self.dom = dom
        if isinstance(cod, dict):
            if cod == {}:
                raise ValueError("need non empty dict argument")
            elif cod[x] == 0 and cod[y] == 0 and cod[z] == 0:
                raise ValueError("(0, 0, 0) is not a P2 point")
            self.cod = {x: cod[x], y: cod[y], z: cod[z]}
        elif isinstance(cod, list):
            if len(cod) == 0:
                raise ValueError("need non empty list")
            elif cod == [0, 0, 0]:
                raise ValueError("(0, 0, 0) is not a P2 point")
            self.cod = {x: cod[0], y: cod[1], z: cod[2]}
        elif isinstance(cod, tuple):
            if len(cod) == 0:
                raise ValueError("need non empty tuple")
            elif cod == (0, 0, 0):
                raise ValueError("(0, 0, 0) is not a P2 point")
            self.cod = {x: cod[0], y: cod[1], z: cod[2]}

    def __repr__(self):
        return 'P2Point([%s: %s: %s], %s)' % (c for c in self.cod.values()) + (self.dom.as_SFF())

    def __str__(self):
        return 'P2Point([%s: %s: %s], %s)' % (c for c in self.cod.values()) + (self.dom.as_SFF())

    def __eq__(p, q):
        raise NotImplementedError

def sff(rep, mod):
    """Constructor method for SFF"""
    return SFF(rep, mod)

def sffpoly(rep, dom):
    """Constructor method for SFFPoly"""
    return SFFPoly(rep, dom)

def p2point(cod, dom):
    """Constructor method for SFFPoly"""
    return P2Point(cod, dom)

def reduce(f, dom):
    if not isinstance(f, Expr):
        raise TypeError("reduce() argument must be an Expr object, not %s", type(f))
    elif isinstance(f, Integer) or isinstance(f, int):
        return f
    var = poly(f).gens
    for rel in dom.rel_list:
        if rel['rep'] == 0 or not rel['var'] in var:
            continue
        f = simple_reduce(f, rel)
    if isinstance(f, Integer) or isinstance(f, int):
        return f % dom.mod
    return poly(f, domain=dom.as_sympy_FF()).as_expr()

def simple_reduce(f, rel):
    if not isinstance(f, Expr):
        raise TypeError("reduce() argument must be an Expr object, not %s", type(f))

    _lm = rel['var'] ** rel['deg']
    if poly(f).degree(rel['var']) < poly(rel['rep']).degree(rel['var']):
        return f
    else:
        return _simple_reduce(f, _lm, _lm - rel['rep'], rel['var'])

def _simple_reduce(f, lm, sub, var):
    if poly(f).degree(var) == poly(lm).degree(var):
        return expand(f.subs({lm: sub}))
    else:
        return expand(_simple_reduce(f, lm * var, sub * var, var).subs({lm: sub}))

def ff_solve(polys):
    if len(polys) == 0:
        raise TypeError("solve() argument must be one or more sffpolys")
    # _vars = list(set(_v for _v in _p.var for _p in polys))
    # _mindeg = min(item.degree() for item in polys)
    _stack = polys[0].simple_solve()
    _sol = _stack
    for p in polys:
        if p == polys[0]:
            continue
        for q in _stack:
            if not reduce(p.subs(q), polys[0].dom) == 0:
                _sol.remove(q)
    return _sol