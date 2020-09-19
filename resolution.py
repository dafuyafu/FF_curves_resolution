'''

	Developer's note (07/31/2020): 
	We will at first implemtent the outline of the Algorithm S-M* without any validation.

'''
from itertools import product

from sympy.core.expr import Expr
from sympy.core.function import diff
from sympy.core.symbol import symbols
from sympy.polys.polytools import Poly, poly, resultant, factor_list

from sffpoly import sff, sffpoly, ff_solve, reduce

def sing(f):
	if not isinstance(f, Poly):
		if isinstance(f, SFFPoly):
			f = poly(f.rep, domain=f.as_sympy_FF())
		else:
			raise TypeError("argument must be a Poly or SFFPoly object, not %s" % f.__class__.__name__)
	x, y = symbols('x y') # tuple
	mod = f.get_modulus()
	p = factor_list(resultant(f, diff(f, x), y), modulus=mod)
	q = factor_list(resultant(f, diff(f, y), x), modulus=mod)
	_list = []

	for _p in p[1]:
		if _has_roots(_p[0], x, mod):
			pass
		
	if not _has_roots(p, x, mod):
		_p = p.subs({x: a})
		_list.append(_p)
		_flag_p = True
	else:
		_p = 0
		_flag_p = False
	if not _has_roots(q, y, mod):
		_q = q.subs({y: b})
		_list.append(_q)
		_flag_q = True
	else:
		_q = 0
		_flag_q = False

	_sff = sff(_list, mod)
	_sff_p = sff(_p, mod)
	_sff_q = sff(_q, mod)
	f = sffpoly(f, _sff)
	p = sffpoly(p, _sff_p)
	q = sffpoly(q, _sff_q)

	if _flag_p:
		_sol_p = [{x: reduce(a ** (p.dom.mod ** e), _sff_p)} for e in range(p.degree())]
	else:
		_sol_p = [{x: n} for n in range(p.dom.mod) if p.subs({x: n}) == 0]
	if _flag_q:
		_sol_q = [{y: reduce(b ** (q.dom.mod ** e), _sff_q)} for e in range(q.degree())]
	else:
		_sol_q = [{y: n} for n in range(q.dom.mod) if q.subs({y: n}) == 0]
	_sol = []

	for point in product(_sol_p, _sol_q):
		_point = {x: point[0][x], y: point[1][y]}
		if p.subs(_point) == 0 and q.subs(_point) == 0 and f.subs(_point) == 0:
			_sol.append(_point)
	return _sol, _sff.as_SFF()

def _has_roots(f, var, mod):
	for i in range(mod):
		if f.subs({var: i}) % mod == 0:
			return True
	return False

def radical(f, mod):
	factors = factor_list(f.as_expr(), modulus=mod)
	radical = 1
	for fac in factors[1]:
		radical *= fac[0]
	return radical