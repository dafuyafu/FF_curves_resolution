from sympy.core.symbol import symbols
from sff import sff
from sffpoly import sffpoly

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