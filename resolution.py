'''

	Developer's note (07/31/2020): 
	We will at first implemtent the outline of the Algorithm S-M* without any validation.

'''

from sympy.core import Expr, symbols
from sympy.polys.polytools import poly, resultant, factor_list
# from sffpoly import SFFPoly, P2Point
from ffpoly import SFF, FFPoly, P2Point

def find_sing(f):
# 	if not isinstance(f, Poly):
# 		if isinstance(f, )

	# validates whether f in k[X,Y,Z]^h

	var = symbols('x,y,z')
	deg = f.degree()
	sing = []
	c = 0
	for v in var:
		_f = f.as_expr()
		for u in var:
			if v == u:
				continue
			_f = _f.subs({u: u * v})
		_f = quo(_f, v ** deg)
		_var = var.remove(v)
		_f = _wrap(_f,f)
		p = [_wrap(resultant(_f, diff(_f, _var[0]), _var[1]), f), 
			 _wrap(resultant(_f, diff(_f, _var[1]), _var[0]), f)]
		sol = []
		'''
			Assume that only one of p[0] and p[1] has nonrational points.
		'''
		for i in range(len(_var)):
			if _has_roots(p[i]):
				sol.append(solve(p[i]))
			else:
				_rel_var = symbols('a_' + str(c))
				c += 1
				_rel = radical(p[i]).as_expr().subs({var[i]: _rel_var})
				sol.append([_rel_var, - _rel_var])
		points = direct_product(sol[0], sol[1])
		



def _has_roots(f):
	if not isinstance(f, Poly):
		raise ValueError("needs to be Poly object")
	if not len(f.gens) == 1:
		raise ValueError("needs to be monovariate polynomial")
	mod = f.get_modulus()
	for i in range(mod):
		if f.subs({f.gens[0]: i}) == 0:
			return True
	return False

def _wrap(ex, dom):
	return poly(ex, domain=dom.domain)

def _wrap_sff(ex, rel, dom):
	return SFFPoly(ex, rel, domain=dom.domain)

def radical(f):
	factors = factor_list(f.as_expr(), modulus=f.get_modulus())
	radical = 1
	for fac in factors[1]:
		radical *= fac[0]
	return _wrap(radical, f)

def direct_product(list1, list2):
	_list = []
	for l1 in list1:
		for l2 in list2:
			_list.append((l1, l2))
	return _list