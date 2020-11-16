'''

	Developer's note (07/31/2020): 
	We will at first implemtent the outline of the Algorithm S-M* without any validation.

'''
from itertools import product
import time

from sympy.core.expr import Expr
from sympy.core.function import diff
from sympy.core.symbol import symbols
from sympy.polys.polytools import Poly, poly, resultant, factor_list

from sffpolytools import sffpoly, reduce
from sffdomains import sff

def sing(f):
	start = time.time()
	if not isinstance(f, Poly):
		if isinstance(f, SFFPoly):
			f = poly(f.rep, domain=f.as_sympy_FF())
		else:
			raise TypeError("argument must be a Poly or SFFPoly object, not %s" % f.__class__.__name__)
	x, y = symbols('x y') # tuple
	mod = f.get_modulus()
	f_x = factor_list(resultant(f, diff(f, x), y), modulus=mod)[1]
	f_y = factor_list(resultant(f, diff(f, y), x), modulus=mod)[1]

	rel_list = []
	count, a, p_x, sol_x = 0, [], [], []
	for f_x_ in f_x:
		f_x_ = f_x_[0]
		if not _has_roots(f_x_, x, mod):
			a.append(symbols('a_' + str(count)))
			f_x_a = f_x_.subs({x: a[count]})
			rel_list.append(f_x_a)
			p_x_ = sffpoly(f_x_, sff(f_x_a, mod))
			p_x.append(p_x_)
			sol_x.extend([{x: reduce(a[count] ** (mod ** d), p_x_.dom)} for d in range(p_x_.degree())])
			count += 1
		else:
			p_x_ = sffpoly(f_x_, sff(0, mod))
			p_x.append(p_x_)
			sol_x.extend([{x: n} for n in range(mod) if p_x_.subs({x: n}) == 0])

	count, b, p_y, sol_y = 0, [], [], []
	for f_y_ in f_y:
		f_y_ = f_y_[0]
		if not _has_roots(f_y_, y, mod):
			b.append(symbols('b_' + str(count)))
			f_y_b = f_y_.subs({y: b[count]})
			rel_list.append(f_y_b)
			p_y_ = sffpoly(f_y_, sff(f_y_b, mod))
			p_y.append(p_y_)
			sol_y.extend([{y: reduce(b[count] ** (mod ** d), p_y_.dom)} for d in range(p_y_.degree())])
			count += 1
		else:
			p_y_ = sffpoly(f_y_, sff(0, mod))
			p_y.append(p_y_)
			sol_y.extend([{y: n} for n in range(mod) if p_y_.subs({y: n}) == 0])

	sff_ = sff(rel_list, mod)
	f = sffpoly(f, sff_)
	f_x = f.diff(x)
	f_y = f.diff(y)
	sol_f = []
	for point in product(sol_x, sol_y):
		point_ = {x: point[0][x], y: point[1][y]}
		if f_x.subs(point_) == 0 and f_y.subs(point_) == 0 and f.subs(point_) == 0:
			sol_f.append(point_)

	elapsed_time = time.time() - start
	print("elapsed_time:{0}".format(elapsed_time) + "[sec]")
	return sol_f, sff_.as_SFF()

def sing_apart(f):
	if not isinstance(f, Poly):
		if isinstance(f, SFFPoly):
			f = poly(f.rep, domain=f.as_sympy_FF())
		else:
			raise TypeError("argument must be a Poly or SFFPoly object, not %s" % f.__class__.__name__)
	x, y, a, b = symbols('x y a b') # tuple
	mod = f.get_modulus()
	f_x = factor_list(diff(f, x), modulus=mod)[1][0][0].as_expr()
	f_y = factor_list(diff(f, y), modulus=mod)[1][0][0].as_expr()

	rel_list = []
	if not _has_roots(f_x, x, mod):
		f_x_ = f_x.subs({x: a})
		rel_list.append(f_x_)
		f_x = sffpoly(f_x, sff(f_x_, mod))
		sol_x = [{x: reduce(a ** (mod ** d), f_x.dom)} for d in range(f_x.degree())]
	else:
		f_x = sffpoly(f_x, sff(0, mod))
		sol_x = [{x: n} for n in range(f_x.dom.mod) if f_x.subs({x: n}) == 0]

	if not _has_roots(f_y, y, mod):
		f_y_ = f_y.subs({y: b})
		rel_list.append(f_y_)
		f_y = sffpoly(f_y, sff(f_y_, mod))
		sol_y = [{y: reduce(b ** (mod ** d), f_y.dom)} for d in range(f_y.degree())]
	else:
		f_y = sffpoly(f_y, sff(0, mod))
		sol_y = [{y: n} for n in range(f_y.dom.mod) if f_y.subs({y: n}) == 0]

	sff_f = sff(rel_list, mod)
	f = sffpoly(f, sff_f)
	sol_f = []
	for point in product(sol_x, sol_y):
		point_ = {x: point[0][x], y: point[1][y]}
		if f_x.subs(point_) == 0 and f_y.subs(point_) == 0 and f.subs(point_) == 0:
			sol_f.append(point_)
	return sol_f, sff_f.as_SFF()

def _has_roots(f, var, mod):
	for i in range(mod):
		if f.subs({var: i}) % mod == 0:
			return True
	return False