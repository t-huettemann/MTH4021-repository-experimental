# [[file:../README.org::*Class: abstract_abelian_group][Class: abstract_abelian_group:1]]
from gmpy2 import mpz

class abstract_abelian_group:
    """This is the base class for abelian groups. This class contains no data, but defines a setup method that should be called by subclasses during initialisation."""
    
    def __init__(self):
        self.element = globals()[self.__class__.__name__ + '_element']

    def el(self, v):
        return self.element(v, self)
        
    def zero(self):
        return self.element(0, self)

    def to_string(self, v):
        return str(v)

    def equals(self, a, b):
        return self.is_zero(self.sub(a,b))

    def __repr__(self):
        return self.__class__.__name__ + '()'

    def __str__(self):
        return 'abstract abelian group'

    def __eq__(self, b):
        return self.__repr__() == b.__repr__()  # overload "=="

    def __ne__(self, b):       # overload "!="
        return self.__repr__() != b.__repr__()
    

class abstract_abelian_group_element:
    """This is the base class for elements of an abelian group. This class defines generic methods for the algebraic structure."""

    def __init__(self, v, G):
        self.group  = G
        self.value = G.normalise(v)

    def __str__(self):
        return self.group.to_string(self.value)

    def __repr__(self):
        return self.__class__.__name__ + '(' + self.value.__repr__() + ', ' + self.group.__repr__() + ')'
        
    def add(self, b):
        if (self.group != b.group):
            print("Addition in abelian groups: Summands must lie in same group")
            raise
        return self.__class__(self.group.add(self.value, b.value), self.group)

    def sub(self, b):
        # need to implement type checking for b
        if (self.group != b.group):
            print("Subtraction in abelian groups: Arguments must lie in same group")
            raise
        return self.__class__(self.group.sub(self.value, b.value), self.group)

    def is_zero(self):
        return self.group.is_zero(self.value)
    
    def equals(self, b):
        return self.group.equals(self.value, b.value)

    def __add__(self, b):      # overload "+"
        return self.add(b)
        
    def __sub__(self, b):      # overload "-"
        return self.sub(b)
    
    def __eq__(self, b):
        return self.__repr__() == b.__repr__()  # overload "=="

    def __ne__(self, b):       # overload "!="
        return self.__repr__() != b.__repr__()

    def __neg__(self):         # overload unary "-"
        return self.group.zero().sub(self)

    def __iadd__(self, b):
        return self.add(b)

    def __isub__(self, b):
        return self.sub(b)

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, b):
        return self.__repr__() == b.__repr__()
    
    def __ne__(self, b):
        return self.repr() != b.__repr__()

    def power(self, k):  # reasonably fast exponentiation
        # self is a group element
        # k is an integer (possibly mpz)
        r = self.group.zero()
        e = mpz(k)
        if e == 0:
            return r
        b = self
        if b.is_zero():
            return r
        if e < 0:
            b = -b
            e = -e
        for i in reversed(e.digits(2)):
            if i == '1':
                r = r + b
            b = b+b
        return r

    def __mul__(self, k):
        return self.power(k)

    def __rmul__(self, k):
        return self.power(k)
# Class: abstract_abelian_group:1 ends here
