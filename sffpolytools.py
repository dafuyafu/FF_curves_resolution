from itertools import product
from multiprocessing import Pool, Queue, Manager
from queue import Empty
import os

from sympy.core.numbers import Integer
from sympy.core.expr import Expr
from sympy.core.function import diff, expand
from sympy.core.symbol import symbols
from sympy.ntheory.primetest import isprime
from sympy.polys.polytools import factor_list, LC, LT, Poly, poly, resultant

from sffdomains import sff, SFF
from multiprocessingtools import SFFPool, StopEval

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
        
        if isinstance(rep, Integer) or isinstance(rep, int):
            self.rep = rep
            self.var = []
            self.is_int = True
        else:
            self.rep = reduce(rep.as_expr(), dom)
            self.var = [v for v in poly(self.rep).gens if not v in dom.var_list] 
            self.is_int = False

        self.dom = dom
        if len(self.var) == 0:
            self.is_const = True
            self.is_uni = False
        elif len(self.var) == 1:
            self.is_const = False
            self.is_uni = True
        else:
            self.is_const = False
            self.is_uni = False

    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__.__name__, self.rep, self.dom.as_SFF())

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
        if not isinstance(g, SFFPoly):
            if isinstance(g, int) or isinstance(g, Integer):
                return sffpoly(reduce(f.rep + g, f.dom), f.dom)
            else:
                raise TypeError("cannot add %s and %s" % (f.__class__.__name__, g.__class__.__name__))
        if f.dom == g.dom:
            _add = reduce(f.rep + g.rep, f.dom)
            return sffpoly(_add, f.dom)
        else:
            raise ValueError("argument sffpolys have different domains")

    def __sub__(f,g):
        """
        subtract two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = sffpoly(2 * x - a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = sffpoly(x + 4 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f - g
        Ssffpoly(x, x, a, modulus=5)
        """
        if not isinstance(f, SFFPoly) or not isinstance(g, SFFPoly):
            if isinstance(g, int) or isinstance(g, Integer):
                return sffpoly(reduce(f.rep - g, f.dom), f.dom)
            else:
                raise TypeError("cannot add %s and %s" % (f.__class__.__name__, g.__class__.__name__))
        if f.dom == g.dom:
            _sub = reduce(f.rep - g.rep, f.dom)
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
        if isinstance(g, int) or isinstance(g, Integer):
        	return sffpoly(f.rep * g, f.dom)
        if isinstance(g, SFFPoly):
	        if f.dom == g.dom:
	            _mul = reduce(f.rep * g.rep, f.dom)
	            return sffpoly(_mul, f.dom)
	        else:
	            raise ValueError("argument sffpolys have different domains")
        else:
        	raise TypeError("cannot multiple %s and %s" % (f.__class__.__name__, g.__class__.__name__))

    def __truediv__(f, g):
        raise TypeError("cannot divide f by g")

    def __floordiv__(f, g):
    	if f.is_const:
    		return g * f ** (f.dom.num - 2)
    	if g.is_const:
    		return f * g ** (f.dom.num - 2)
    	if not f.is_uni or not g.is_uni:
    		raise TypeError("cannot divide multivariate polynomial(s)")
    	if not f.var[0] == g.var[0]:
    		raise ValueError("cannot divide polynomials which have different variables")
    	if f.degree() < g.degree():
    		return sffgen(0, f.dom)
    	_var = f.var[0]
    	_div = sffpoly(0, f.dom)
    	while f.degree() >= g.degree():
    		lt_f = sffconst(lc(f, _var), f.dom)
    		lt_g = sffconst(lc(g, _var), f.dom)
    		_term = (lt_f / lt_g).rep * _var ** (f.degree() - g.degree())
    		_div += sffpoly(_term, f.dom)
    		_mul = reduce(_term * g.rep, f.dom)
    		f -= sffpoly(_mul, f.dom)
    		if f == 0:
    			break
    	return _div

    def __mod__(f, g):
    	return f - (f // g) * g

    def __pow__(f, e):
        if not isinstance(e, int) and not isinstance(e, Integer):
            raise TypeError("second argument needs to be an integer, not %s" % e.__class__.__name__)
        if f.is_int:
            return sffpoly(f.rep ** e, f.dom)
        if e < 0:
            e = f.dom.num + e - 1
        if e == 0:
            return sffpoly(1, f.dom)
        if e == 1:
            return f
        num_ = bin(e).replace('0b','')
        len_ = len(num_)
        list_ = [(f, len_ - d - 1) for d in range(len_) if num_[d] == '1']
        with Pool(os.cpu_count()) as p:
            result = p.starmap(_pow_self, list_)
        pow_ = 1
        for r in result:
            pow_ = reduce(pow_ * r, f.dom)
        return sffpoly(pow_, f.dom)

    def __eq__(f,g):
        if isinstance(g, SFFPoly):
            if f.rep == g.rep:
                return True
            else:
                return False
        else:
            if f.rep == g:
                return True
            else:
                return False

    def as_expr(self):
        return self.rep

    def as_poly(self):
        if self.is_int:
            raise ValueError("this is an integer.")
        else:
            return poly(self.rep, domain=self.dom.as_sympy_FF())

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

    def get_modulus(self):
        return self.dom.mod

    def subs(self, point):
        """
        substitute points to its variables
        """
        if self.is_int:
            raise TypeError("this is a modular integer")
        if point == {}:
            raise ValueError
        _rep = self.rep
        for k, v in point.items():
            _rep = expand(_rep.subs({k: v}))
        return reduce(_rep, self.dom)

    def subs_as_sffpoly(self, **args):
        return sffpoly(self.subs(args), self.dom)
        
    def solve_abs(self):
        """ solve self over its algebraic closure """
        raise NotImplementedError

    def solve(self):
        # _sol = []
        # for point in self.dom.points_as_dict_iter(self.var):
        #     if self.subs(point) == 0:
        #         _sol.append(point)
        # return _sol
        _sol = []
        m = Manager()
        q = m.Queue()
        with Pool(os.cpu_count()) as p:
        	p.starmap_async(_eval_solve, self.dom.points_as_dict_with_poly_and_queue_iter(self.var, self, q))
        	if q.qsize == self.degree():
        		p.terminate()
        while not q.empty():
        	_sol.append(q.get())
        return _sol

    def is_primitive():
        raise TypeError("This is SFFPoly object.")

    def is_const(self):
        return self.is_const

    def diff(self, var):
        return sffpoly(reduce(diff(self.rep, var), self.dom), self.dom)

    def sing(self):
        """ find singular locus of self """
        pass

    def reduce(self):
        self.rep = reduce(self.rep, self.dom)

    def diff(self, *gens):
    	if len(gens) == 0:
    		return sffpoly(diff(self.rep, self.var[0]), self.dom)
    	else:
    		return sffpoly(diff(self.rep, gens[0]), self.dom)

    def toSFFConst(self):
    	if self.is_const:
    		return sffconst(self.rep, self.dom)
    	else:
    		raise TypeError("cannot convert not constant elements to SFFConst")

    @classmethod
    def random_poly(cls, var, dom, deg, quo=0):
        p = 0
        for i in range(deg):
            p += dom.rand() * var ** i
        if quo == 0:
            return cls(p, dom)
        else:
            return cls(p, dom, quo)

class SFFQuotientPoly(SFFPoly):
    def __init__(self, rel, dom, quo): 
        super().__init__(rel, dom)
        mod = dom.mod
        if quo == 0:
            pass
        else:
            self.quo_list = []
            if isinstance(quo, Expr):
                if not LC(quo.as_poly()) == 1:
                    quo = poly(quo.as_expr() * pow(LC(quo.as_poly()), mod - 2, mod), domain=dom.as_sympy_FF()).as_expr()
                if len(poly(quo).gens) == 1:
                    self.quo_list.append({'var': poly(quo).gens[0],
                                          'rep': quo.as_expr(), 
                                          'deg': poly(quo).degree()})
                else:
                    self.quo_list.append({'var': poly(quo).gens[0], 
                                          'rep': quo.as_expr(), 
                                          'deg': poly(quo).degree()})
            elif isinstance(quo, list):
                if isinstance(quo[0], dict):
                    self.quo_list = quo
                else:
                    for _p in quo:
                        if _p == 0:
                            continue
                        if not LC(_p.as_poly()) == 1:
                            _p = poly(_p * pow(LC(_p.as_poly()), mod - 2, mod), domain=dom.as_sympy_FF()).as_expr()
                        if len(poly(_p).gens) == 1:
                            self.quo_list.append({'var': poly(_p).gens[0], 
                                                  'rep': _p.as_expr(), 
                                                  'deg': poly(_p).degree()})
                        else:
                            self.quo_list.append({'var': poly(_p).gens[0], 
                                                  'rep': _p.as_expr(), 
                                                  'deg': poly(_p).degree()})
            else:
                raise ValueError("the third argument needs to be 0, an Expr instance or list.")

    def __repr__(self):
        return "%s(%s, %s, over %s)" % (self.__class__.__name__, self.rep, self.dom.as_SFF(), self.quos())

    def __add__(f, g):
        add = super().__add__(g)
        return sffquotientpoly(add.as_expr(), f.dom, f.quo_list)

    def __sub__(f, g):
        sub = super().__sub__(g)
        return sffquotientpoly(sub.as_expr(), f.dom, f.quo_list)

    def __mul__(f, g):
        """
        multiple two polynomials ``f`` and ``g``

        Examples
        ========

        >>> f = Ssffpoly(x - a * y, a ** 2 - 2, domain='FF(5)')
        >>> g = Ssffpoly(x - 2 * a * y, a ** 2 - 2, domain='FF(5)')
        >>> f * g
        Ssffpoly(x ** 2 + 2 * a * x * y - y ** 2, x, modulus=5)
        """
        if isinstance(g, int) or isinstance(g, Integer):
            return sffpoly(f.rep * g, f.dom)
        elif isinstance(g, SFFQuatientPoly):
            if f.dom == g.dom:
                _mul = reduce(f.rep * g.rep, f.dom)
                for q in f.quo_list:
                    _mul = simple_reduce(_mul, q)
                return sffpoly(_mul, f.dom)
            else:
                raise ValueError("argument sffpolys have different domains")
        elif isinstance(g, SFFPoly):
            if f.dom == g.dom:
                _mul = reduce(f.rep * g.rep, f.dom)
                return sffpoly(_mul, f.dom)
            else:
                raise ValueError("argument sffpolys have different domains")
        else:
            raise TypeError("cannot multiple %s and %s" % (f.__class__.__name__, g.__class__.__name__))

    def __pow__(f, e):
        if not isinstance(e, int) and not isinstance(e, Integer):
            raise TypeError("second argument needs to be an integer, not %s" % e.__class__.__name__)
        if f.is_int:
            return sffquotientpoly(f.rep ** e, f.dom, f.quo_list)
        if e < 0:
            e = f.dom.num + e - 1
        if e == 0:
            return sffquotientpoly(1, f.dom, f.quo_list)
        if e == 1:
            return f
        num_ = bin(e).replace('0b','')
        len_ = len(num_)
        list_ = [(f, len_ - d - 1) for d in range(len_) if num_[d] == '1']
        with Pool(os.cpu_count()) as p:
            result = p.starmap(_pow_self_quo, list_)
        pow_ = 1
        for r in result:
            pow_ = reduce(pow_ * r, f.dom)
            for q in f.quo_list:
                pow_ = simple_reduce(pow_, q)
        return sffquotientpoly(reduce(pow_, f.dom), f.dom, f.quo_list)

    def quos(self):
        return [q['rep'] for q in self.quo_list]

class SFFConst(SFFPoly):
    def __truediv__(f,g):
        return f * g ** (f.dom.num - 2)

    def is_primitive(self):
        num_ = (self.dom.num - 1) // 2
        if not self ** num_ == -1:
            return False
        with SFFPool(os.cpu_count()) as p:
            try:
                p.starmap_async(_is_primitive, self.iter_with_factor_count(num_)).get()
            except StopEval:
                p.close()
                return False
        return True

    def iter_with_count(self, m):
        for i in range(m):
            yield (self, i)

    def iter_with_factor_count(self, m):
        for i in range(1, m):
            if m % i == 0:
                yield (self, i)
            else:
                continue

    def minpoly(self, var):
        p = sffpoly(1, self.dom)
        var = sffpoly(var, self.dom)
        for i in range(self.dom.exdeg):
            p *= (var - self ** (self.dom.mod ** i))
        return reduce(p.rep, self.dom)

    def toSFFConst(self):
    	raise TypeError("this is already SFFConst")

class SFFInt(SFFConst):
    def is_primitive(self):
        if self.dom.is_prime:
            if not self.rep == 0 and not self.rep == 1:
                return True
            else:
                return False
        else:
            return False

    def toSFFConst(self):
    	raise TypeError("this is already SFFConst")

def sffgen(rep, dom, quo=0):
    if quo == 0:
        if isinstance(rep, int) or isinstance(rep, Integer):
            return sffint(rep, dom)
        elif isintance(rep, Expr):
            var = [v for v in rep.as_poly().gens if not v in dom.var_list]
            if len(var) == 0:
                return sffconst(rep, dom)
            else:
                return sffpoly(rep, dom)
    else:
        return sffquotientpoly(rep, dom, quo)

def sffpoly(rep, dom):
    """Constructor method for SFFPoly"""
    return SFFPoly(rep, dom)

def sffconst(rep, dom):
    return SFFConst(rep, dom)

def sffint(rep, dom):
    return SFFInt(rep, dom)

def sffquotientpoly(rep, dom, quo):
    return SFFQuotientPoly(rep, dom, quo)

def reduce(f, dom):
	if not isinstance(f, Expr):
		raise TypeError("reduce() argument must be an integer or an Expr object, not %s" % f.__class__.__name__)
	elif isinstance(f, Integer) or isinstance(f, int):
		f %= dom.mod
		if f > dom.mod // 2:
			return f - dom.mod
		else:
			return f
	else:
		var = poly(f).gens
		f = expand(f)
		for rel in dom.rel_list:
			if rel['rep'] == 0 or not rel['var'] in var:
				continue
			f = simple_reduce(f, rel)
		if isinstance(f, Integer) or isinstance(f, int):
			f %= dom.mod
			if f > dom.mod // 2:
				return f - dom.mod
			else:
				return f
		else:
			return poly(f, domain=dom.as_sympy_FF()).as_expr()

def simple_reduce(f, rel):
    if not isinstance(f, Expr):
        raise TypeError("reduce() argument must be an integer or Expr object, not %s" % f.__class__.__name__)
    lm = rel['var'] ** rel['deg']
    if poly(f).degree(rel['var']) < poly(rel['rep']).degree(rel['var']):
        return expand(f)
    else:
        return _simple_reduce(f, lm, lm - rel['rep'], rel['var'])

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

def _eval_solve(f, p, q):
	if f.subs(p) == 0:
		q.put(p)

def _pow_self(f, n):
    """ f needs to be an sffpoly element"""
    pow_ = f.rep
    for i in range(n):
        pow_ = reduce(pow_ * pow_, f.dom)
    if isinstance(pow_, int) or isinstance(pow_, Integer):
        return pow_
    else:
        return poly(pow_, domain=f.dom.as_sympy_FF()).as_expr()

def _pow_self_quo(f, n):
    """ f needs to be an sffpoly element"""
    pow_ = f.rep
    for i in range(n):
        pow_ = reduce(pow_ * pow_, f.dom)
        for q in f.quo_list:
            pow_ = simple_reduce(pow_, q)
    if isinstance(pow_, int) or isinstance(pow_, Integer):
        return pow_
    else:
        return poly(pow_, domain=f.dom.as_sympy_FF()).as_expr()

def _is_primitive(f, i):
    if f ** i == -1:
        raise StopEval()
    else:
        pass

def primitive_elements(dom):
    result = []
    for e in dom.elements_iter():
        if sffconst(e, dom).is_primitive():
            result.append(e)
    return result

def primitive_element(dom):
    for e in dom.elements_iter():
        if sffconst(e, dom).is_primitive():
            return sffconst(e, dom)
    raise TypeError("this has no primitive element")

def _primitive_eval(f, q):
    if f.is_primitive():
        q.put(f.rep)

def simplify(ff, var):
    if not isinstance(ff, SFF):
        raise TypeError("argument must be a SFF object. not %s" % ff.__class__.__name__)
    minpoly_ = primitive_element(ff).minpoly(var)
    return sff(minpoly_, ff.mod)

def lc(f, *gens):
	if len(gens) == 0:
		return LC(poly(f.rep), f.var[0])
	else:
		return LC(poly(f.rep), gens[0])