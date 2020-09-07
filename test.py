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
	x, y, a, b = symbols('x, y, a, b') # tuple
	mod = f.get_modulus()
	p = radical(resultant(f, diff(f, x), y), mod)
	q = radical(resultant(f, diff(f, y), x), mod)
	_list = []

	if not _has_roots(p, x, mod):
		_p = p.subs({x: a})
		_list.append(_p)
	else:
		_p = 0
	if not _has_roots(q, y, mod):
		_q = q.subs({y: b})
		_list.append(_q)
	else:
		_q = 0

	_sff = sff(_list, mod)
	_sff_p = sff(_p, mod)
	_sff_q = sff(_q, mod)
	_f = sffpoly(f, _sff)
	_p = sffpoly(p, _sff_p)
	_q = sffpoly(q, _sff_q)
	_sol_p = _p.simple_solve()
	_sol_q = _q.simple_solve()
	_sol = []

	for point in product(_sol_p, _sol_q):
		_point = {x: point[0][x], y: point[1][y]}
		if reduce(_p.subs(_point), _sff) == 0 and reduce(_q.subs(_point), _sff) == 0:
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