from sympy.polys.polytools import Poly
from sffpolytools import sffpoly, lc, SFFPoly, SFFQuotientPoly, sffquotientpoly

def sffgcd(f, g):
	q = f % g
	if q == 0:
		return g
	else:
		return sffgcd(g, q)

def sffsff_list(f):
	if not f.is_uni:
		raise TypeError("cannot sffsff multivariate polynomial")
	if f == sffpoly(0, f.dom):
		return f
	_lc = lc(f) ** (f.dom.mod - 2) % f.dom.mod
	f *= _lc
	fact = sffpoly(1, f.dom)
	fact_list = []
	flat = f // sffgcd(f, f.diff())
	m = 0
	while not flat.is_const:
		while f % flat == 0:
			f = f // flat
			m += 1
		_flat = sffgcd(flat, f)
		g = flat // _flat
		flat = _flat
		fact *= g ** m
		fact_list.append((g, m))
	if not f.is_const:
		p = f.dom.mod
		colist = f.as_poly().as_list()
		newcolist = []
		for i in len(colist) / p:
			newcolist.append(colist[i * p])
		newsffpoly = sffpoly(Poly.from_list(newcolist, gens=f.var[0], domain=f.dom.as_sympy_FF()), f.dom)
		_fact_list = sffsff_list(newsffpoly)[1]
		_list = []
		for t in _fact_list:
			_list.append((t[0], t[1] * f.dom.mod))
		fact_list.extend(_list)
	return (_lc, fact_list)

def sfffactor_list(f):
	_sffsff_list = sffsff_list(f)
	_list = []

def sfffactor_list_berlekamp(f):
	if not isinstance(f, SFFPoly):
		raise TypeError("needed a SFFPoly object, not %s" % f.__class__.__name__)
	
def sfffactor_equal_degree(f, d):
	if not isinstance(f, SFFPoly):
		raise TypeError("needed a SFFPoly object, not %s" % f.__class__.__name__)
	r = f.degree() // d
	q = f.dom.num
	F = [f]
	var = f.var[0]
	while len(F) < r:
		g = SFFQuotientPoly.random_poly(var, f.dom, r * d, f.rep)
		g = sffpoly((g ** ((q ** d - 1) // 2) - 1).rep, f.dom)
		F_1 = []
		while len(F) > 0:
			h = F.pop(-1)
			z = sffgcd(h, g)
			if z == 1 or z == h:
				F_1.append(h)
			else:
				F_1.append(z)
				F_1.append(h // z)
		F = F_1
	return F