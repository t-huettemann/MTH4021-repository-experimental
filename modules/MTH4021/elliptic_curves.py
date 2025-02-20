# [[file:../README.org::*Class elliptic_curve][Class elliptic_curve:1]]
import MTH4021.abelian_groups as ab

class elliptic_curve(ab.abstract_abelian_group):
    """An elliptic curve."""

    def __init__(self, p):
        # the elliptic curve define by polynomial p, which must be an object
        # of type polynomialring_over_field_element
        # of the form  x^3 + ax + b
        # NO CHECKING currently done!
        self.curve = p
        self.a = p.coeff(1)
        self.field = p.ring.basering
        self.element = globals()[self.__class__.__name__ + '_element']
        self.points_cache = set()

    def zero(self):
        return self.element((), self)

    def is_zero(self, v):
        return v == ()

    def add(self, x, y): # x and y are pairs of elements; neutral element represented by empty pair ()
        if x == ():
            return y
        elif y == ():
            return x
        elif x[0] == y[0] and x[1] == -y[1]:
            return ()
        elif x[0] != y[0]:
            xdp = ((y[1]-x[1])/(y[0]-x[0])) * ((y[1]-x[1])/(y[0]-x[0])) - x[0] - y[0];
            return ( xdp, ((y[1]-x[1])/(y[0]-x[0]))*(x[0]-xdp) - x[1] )
        else:
            hh = x[0]*x[0]
            h   = (hh + hh + hh + self.a) / (x[1] + x[1])
            xdp = h*h - x[0] - x[0]
            return (xdp, h * (x[0] - xdp) - x[1])

    def sub(self, x, y):
        if y == ():
            return x
        return self.add(x, (y[0], -y[1]))

    def normalise(self, v):
        return v

    def __repr__(self):
        return self.__class__.__name__ + '(' + self.curve.__repr__() + ')'

    def __str__(self):
        return 'elliptic curve given by ' + self.curve.__str__() + ' over ' + self.curve.ring.basering.__str__()

    def points(self):
        if self.points_cache == set():
            squares = {y*y: y for y in self.curve.ring.basering}
            for x in self.curve.ring.basering:
                yq = self.curve.eval(x)
                if yq in squares:
                    self.points_cache.add(self.el((x,  squares[yq])))
                    self.points_cache.add(self.el((x, -squares[yq])))
            self.points_cache.add(self.el('∞'))
        return self.points_cache


class elliptic_curve_element(ab.abstract_abelian_group_element):
    """The base class for elements of an elliptic curve."""

    def __str__(self):
        if self.is_zero():
            return '∞'
        else:
            return '(' + str(self.value[0]) + ', ' + str(self.value[1]) + ')'

    def __init__(self, v, C):
        if v == '∞':
            super().__init__((), C)
        else:
            super().__init__(v, C)
# Class elliptic_curve:1 ends here
