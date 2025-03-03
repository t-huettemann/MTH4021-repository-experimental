# [[file:../README.org::*Rings and fields][Rings and fields:1]]
import gmpy2
from secrets import randbelow
from gmpy2 import mpz, is_zero, gcdext, gcd, divm, is_prime, f_mod
import itertools
import primefac as pf

class ring:
    """This is the base class for a (commutative) ring. This class contains no data, but defines a setup method that should be called by subclasses during initialisation."""

    def __init__(self):
        self.element = globals()[self.__class__.__name__ + '_element']

    def zero(self):
        return self.element(0, self)

    def one(self):
        return self.element(1, self)

    def to_string(self, v):
        return str(v)

    def equals(self, a, b):
        return self.is_zero(self.sub(a,b))

    def power(self, x, k):  # reasonably fast exponentiation
                            # x a ring element
                            # k is an integer (possibly mpz)
        r = self.one().value
        b = x
        e = mpz(k)
        if e == 0:
            return r
        if e < 0:
            print("Cannot invert elements in generic ring")
            return
        if self.is_zero(x):
            return x
        for i in reversed(e.digits(2)):
            #print("\n\nDigit:", i)
            if i == '1':
                r = self.mult(r, b)
            b = self.mult(b, b)
        return r

    def __repr__(self):
        return self.__class__.__name__ + '()'

    def __str__(self):
        return 'ring'

    def __eq__(self, b):
        return self.__repr__() == b.__repr__()  # overload "=="

    def __ne__(self, b):       # overload "!="
        return self.__repr__() != b.__repr__()

    def __call__(self, *args):
        if isinstance(args[0], self.element):
            return args[0]
        return self.element(*args, self)

    def __pow__ (self, x, k):
        return self.power(x, k)


class ring_element:
    """This is the base class for elements of a ring. This class defines generic methods for the algebraic structure."""

    def __init__(self, v, r):
        self.ring  = r
        self.value = r.normalise(v)

    def __str__(self):
        return self.ring.to_string(self.value)

    def __repr__(self):
        return self.__class__.__name__ + '(' + self.value.__repr__() + ', ' + self.ring.__repr__() + ')'

    def add(self, b):
        if (self.ring != b.ring):
            print("Addition: Summands must lie in same ring")
            raise
        return self.__class__(self.ring.add(self.value, b.value), self.ring)

    def sub(self, b):
        # need to implement type checking for b
        if (self.ring != b.ring):
            print("Subtraction: Arguments must lie in same ring")
            raise
        return self.__class__(self.ring.sub(self.value, b.value), self.ring)

    def mult(self, b):
        # need to implement type checking for b
        if (self.ring != b.ring):
            print("Multiplication: Factors must lie in same ring")
            raise
        return self.__class__(self.ring.mult(self.value, b.value), self.ring)

    def is_zero(self):
        return self.ring.is_zero(self.value)

    def is_one(self):
        return self.__repr__() == self.ring.one().__repr__()

    def equals(self, b):
        return self.ring.equals(self.value, b.value)

    def power(self, k):
        return self.ring(self.ring.power(self.value, k))

    def __add__(self, b):      # overload "+"
        return self.add(b)

    def __sub__(self, b):      # overload "-"
        return self.sub(b)

    def __mul__(self, b):      # overload "*"
        return self.mult(self.ring(b))

    def __rmul__(self, b):      # overload "*"
        return self.mult(self.ring(b))

    def __eq__(self, b):
        return self.__repr__() == b.__repr__()  # overload "=="

    def __ne__(self, b):       # overload "!="
        return self.__repr__() != b.__repr__()

    def __neg__(self):         # overload unary "-"
        return self.ring.zero().sub(self)

    def __iadd__(self, b):
        return self.add(b)

    def __isub__(self, b):
        return self.sub(b)

    def __imul__(self, b):
        return self.mult(b)

    def __pow__(self, k):
        return self.ring(self.ring.power(self.value, k))

    def __hash__(self):
        return hash(self.__repr__())

    def _print_sign(self):
        return '+'

    def _print_abs(self):
        return str(self)

    
class Z(ring):
    def __init__(self):
        # self.modulus = mpz(0)
        # self.iterator = globals()['zmod_iterator']
        self.char = mpz(0) # self.modulus
        self.cardinality = float('inf') # self.modulus
        super().__init__()

    def __str__(self):
        return 'Z'

    def __repr__(self):
        return self.__class__.__name__ + '()' #'(' + self.modulus.__repr__() + ')'

    def normalise(self, v):
        return mpz(v)

    def add(self, a, b):
        return mpz(a) + mpz(b)

    def sub(self, a, b):
        return mpz(a)-mpz(b)

    def mult(self, a, b):
        return mpz(a)*mpz(b)

    def is_zero(self, v):
        return mpz(v) == 0

    def mod(self, m, n):
        return gmpy2.f_mod(mpz(m), mpz(n))
    
    def div(self, m, n):
        return gmpy2.f_div(mpz(m), mpz(n))
    
    # def __iter__(self):
    #     return self.iterator(self)

    # def random_element(self):
    #     """Return a random element."""
    #     return self(randbelow(self.modulus))


class Z_element(ring_element):
    # definition of an element in Z
    # use generic methods from ring_element!

    def _print_sign(self):
        if self.value < 0:
            return '-'
        else:
            return '+'

    def _print_abs(self):
        return str(abs(self.value))

    def mod(self, n):
        return self.ring(self.ring.mod(self.value, n.value))

    def div(self, n):
        return self.ring(self.ring.div(self.value, n.value))
class zmod(ring):
    def __init__(self, m):
        self.modulus = abs(mpz(m))
        self.iterator = globals()['zmod_iterator']
        self.char = self.modulus
        self.cardinality = self.modulus
        super().__init__()

    def __str__(self):
        return 'Z/' + str(self.modulus)

    def __repr__(self):
        return self.__class__.__name__ + '(' + self.modulus.__repr__() + ')'

    def normalise(self, v):
        if isinstance(v, zmod_element):
            return v
        return f_mod(mpz(v), self.modulus)

    def add(self, a, b):
        return f_mod(a + b, self.modulus)

    def sub(self, a, b):
        return f_mod(a - b, self.modulus)

    def mult(self, a, b):
        return f_mod(a * b, self.modulus)

    def is_zero(self, v):
        return (f_mod(mpz(v), self.modulus) == 0)

    def __iter__(self):
        return self.iterator(self)

    def random_element(self):
        """Return a random element."""
        return self(randbelow(self.modulus))


class zmod_element(ring_element):
    # definition of an element in Z/n
    # use generic methods from ring_element!
    pass


class zmod_iterator:
    def __init__(self, ring):
        self.ring = ring
        self.current = self.ring.zero()
        self.done    = False

    def __next__(self):
        if self.done:
            raise StopIteration
        r = self.current
        self.current = r + self.ring.one()
        self.done = self.current.is_zero()
        return r

class field(ring):
    """Base class for fields. Does nothing, and contains no data."""
    def inv(self, x): # compute 1/x
        return self.div(self.one().value, x)

    def power(self, x, k):  # reasonably fast exponentiation
        # x a ring element
        # k is an integer (possibly mpz)
        r = self.one().value
        e = mpz(k)
        if e == 0:
            return r
        b = x
        if e < 0:
            b = self.inv(x)
            e = -e
        if self.is_zero(x):
            return x
        for i in reversed(e.digits(2)):
            if i == '1':
                r = self.normalise(self.mult(r, b))
            b = self.mult(b, b)
        return r


class field_element(ring_element):
    """Base class for elements of a field. Provides method "div" which must be provided by actual implementation of field."""

    def div(self, b):
        if (self.ring != b.ring):
            print("Division: Arguments must lie in same field")
            raise
        if b.is_zero():
            print(f"Method 'div': Cannot divide by 0 (Class {self.__class__})")
        return self.__class__(self.ring.div(self.value, b.value), self.ring)

    def inv(self):
        if self.is_zero():
            print(f"Method 'inv': Cannot invert 0 (Class {self.__class__})")
            raise
        return self.ring.one().div(self)

    def __truediv__(self, b):      # overload "/"
        return self.div(b)

    def __rtruediv__(self, b):
        return self.ring(b).div(self)
    
    def __idiv__(self, b):         # overload "/="
        return self.div(b)

class primefield(zmod, field):

    def __init__(self, p):
        m = abs(mpz(p))
        if is_prime(m):
            self.char = m
            super().__init__(p)
        else:
            print("Modulus of prime field must be prime")
            raise

    def __str__(self):
        return 'F_' + str(self.modulus)

    def div(self, a, b):
        return divm(a, b, self.modulus)


class primefield_element(zmod_element, field_element):
    # definition of a "finite prime field element"
    # should implement checking that modulus is prime
    def __idiv__(self, b):
        return self.div(b)

class Q(field):

    def __init__(self):
        self.char = mpz(0)
        self.cardinality = float('inf')
        self.value = gmpy2.mpq(0)
        super().__init__()

    def __str__(self):
        return 'Q'

    def div(self, a, b):
        return(a/b)

    def mult(self, a, b):
        return(a*b)

    def add(self, a, b):
        return(a+b)

    def sub(self, a, b):
        return(a-b)

    def normalise(self, a):
        return gmpy2.mpq(a)

    def is_zero(self, a):
        return(a == 0)
    
class Q_element(field_element):
    pass

class polynomialring(ring):

    def __init__(self, r, i='x', parentheses=['', '']):
        if isinstance(r, ring):
            self.basering = r
            self.indeterminate = i
            self.parentheses = parentheses
            super().__init__()
        else:
            print("Coefficients must lie in a ring")
            raise

    def __str__(self):
        return f"{self.basering}[{self.indeterminate}]"

    def __repr__(self):
        return self.__class__.__name__ + '(' + self.basering.__repr__() + ', i=\'' + self.indeterminate + '\', parentheses=[\'' + self.parentheses[0] + '\', \'' + self.parentheses[1] + '\'])'

    def normalise(self, v):
        # v is a list with coefficients, starting with degree 0; if it's a single element, turn into a list first
        if isinstance(v, list):
            w = list(map(lambda x: x if isinstance(x, ring_element) else self.basering.element(x, self.basering), v))
        else:
            if isinstance(v, ring_element):
                w = [v]
            else:
                w = [self.basering.element(v, self.basering)]
        while ((w != []) and (w[-1].is_zero())):
            w.pop()
        return w

    def add(self, a, b):
        v = []
        for (s,t) in itertools.zip_longest(a, b, fillvalue=self.basering.zero()):
            v.append(s.add(t))
        return self.normalise(v)

    def sub(self, a, b):
        v = []
        for (s,t) in itertools.zip_longest(a, b, fillvalue=self.basering.zero()):
            v.append(s.sub(t))
        return self.normalise(v)

    def mult(self, a, b):
        da = len(a)-1
        db = len(b)-1
        if da<0:
            return []
        if db<0:
            return []
        r = []
        for d in range(0, da+db+1, 1):
            h = self.basering.zero()
            for i in range(0, min(d, da)+1, 1):
                if d-i > db:
                    pass
                else:
                    h = h.add(a[i].mult(b[d-i]))
            r.append(h)
        return self.normalise(r)

    def zero(self):
        return self.element([], self)

    def one(self):
        return self.element([self.basering.one()], self)

    def to_string(self, v):
        po = self.parentheses[0]
        pc = self.parentheses[1]
        r = []
        sgn = []
        d = len(v)-1
        if d < 0:
            r.append('0')
            sgn.append(' + ')
        else:
            for i in range(d, -1, -1):
                if v[i].is_zero():
                    pass
                elif i == 0:
                    r.append(po + v[0]._print_abs() + pc)
                    #r.append(v[0]._print_abs())
                    sgn.append(' ' + v[0]._print_sign() + ' ')
                elif i == 1:
                    s = v[1]._print_abs()
                    sgn.append(' ' + v[1]._print_sign() + ' ')
                    if s == '1':
                        r.append(self.indeterminate)
                    else:
                        r.append(po + s + pc + self.indeterminate)
                else:
                    s = v[i]._print_abs()
                    sgn.append(' ' + v[i]._print_sign() + ' ')
                    if s == '1':
                        r.append(self.indeterminate + '^' + str(i))
                    else:
                        r.append(po + s + pc + self.indeterminate + '^' + str(i))
        rr = [x for y in zip(sgn, r) for x in y]
        if rr[0] == " + ":
            rr.pop(0)
        elif rr[0] == " - ":
            rr[0] = "-"
        return "".join(rr)
        #     for i in range(d, -1, -1):
        #         if i == 0:
        #             if v[0].is_zero():
        #                 pass
        #             else:
        #                 r.append(po + v[0].__str__() + pc)
        #         elif i == 1:
        #             s = v[1].__str__()
        #             if s == '0':
        #                 pass
        #             elif s == '1':
        #                 r.append(self.indeterminate)
        #             else:
        #                 r.append(po + s + pc + self.indeterminate)
        #         else:
        #             s = v[i].__str__()
        #             if s == '0':
        #                 pass
        #             elif s == '1':
        #                 r.append(self.indeterminate + '^' + str(i))
        #             else:
        #                 r.append(po + s + pc + self.indeterminate + '^' + str(i))
        # return " + ".join(r)

    def deg(self, p):
        q = self.normalise(p)
        return len(q)-1

    def lc(self, p):
        return self.normalise(p)[-1]

    def coeff(self, p, n=0):
        if n < 0:
            return self.basering.zero()
        q = self.normalise(p)
        if n >= len(q):
             return self.basering.zero()
        else:
            return q[n]

    def is_zero(self, p):
        return(self.normalise(p) == [])

    def is_monic(self, p):
        return self.lc(p).is_one()

    def div_mod(self, f,g): # compute q, r with f = gq+r and deg(r) < deg(g). Returns list [q, r].
        if not(self.is_monic(g)):
            print("Class polynomialring: Divisor must me monic for division with remainder")
            raise
        d = self.deg(g)
        r = self.normalise(f)
        q = []
        while self.deg(r) >= d:
            e = self.deg(r) - d
            m = []
            for i in range(0, e, 1):
                m.append(self.basering.zero())
            m.append(self.lc(r))
            q = self.add(q, m)
            r = self.sub(r, self.mult(m, g))
        return [self.element(q, self), self.element(r, self)]

    def div(self, f, g):
        return self.div_mod(f,g)[0]

    def mod(self, f, g):
        return self.div_mod(f,g)[1]
    
    
class polynomialring_element(ring_element):
    # using generic methods mostly!

    def __init__(self, v, r):
        self.ring = r
        self.value = r.normalise(v)

    def deg(self):
        return self.ring.deg(self.value)

    def lc(self):
        return self.ring.lc(self.value)

    def coeff(self, n):
        return self.ring.coeff(self.value, n)

    def eval(self, x):
        if self.ring.basering.__class__.__name__ != x.ring.__class__.__name__:
            print("Class polynomialring_element: To evaluate polynomial, value of variable must lie in coefficient ring")
            raise
        if self.is_zero():
            # return x.__class__.ring.zero()
            return self.ring.basering.zero()
        d = self.deg()
        r = self.coeff(d)
        while d>0:
            d = d-1
            r = r*x + self.coeff(d)
        return r

    def is_monic(self):
        return self.lc().is_one()

    def div(self, g):
        return self.ring.div_mod(self.value, g.value)[0]

    def mod(self, g):
        return self.ring.div_mod(self.value, g.value)[1]

    def __mod__(self, g): # overloading "%"
        return self.mod(g)




class polynomialring_over_field(polynomialring):

    def __init__(self, r, i='x', parentheses=['', '']):
        if isinstance(r,  field):
            super().__init__(r, i, parentheses)
        else:
            print("Class polynomialring_over_field: Coefficients must lie in a field")
            raise

    def div_mod(self, f,g): # compute q, r with f = gq+r and deg(r) < deg(g). Returns list [q, r].
        if self.is_zero(g):
            print("Class polynomialring_over_field: Cannot divide by zero polynomial")
            raise
        d = self.deg(g)
        c = self.lc(g)
        r = self.normalise(f)
        q = []
        while self.deg(r) >= d:
            e = self.deg(r) - d
            m = []
            for i in range(0, e, 1):
                m.append(self.basering.zero())
            m.append(self.lc(r).div(c))
            # print("m = ", str(m))
            # print("q = ", q)
            q = self.add(q, m)
            r = self.sub(r, self.mult(m, g))
        return [q, r]

    def div(self, f, g):
        return self.div_mod(f,g)[0]

    def mod(self, f, g):
        return self.div_mod(f,g)[1]

    def Bezout(self, f,g): # compute a, b so that af + bg = gcd(f,g), with gcd monic (leading coefficient one).
        # Return list [gcd(f,g), a, b].
        # Note: f, g are values (not objects), and the returned list consists of values (not objects).
        if self.is_zero(f):
            return [g, self.zero().value, self.one().value] # changed
        elif self.is_zero(g):
            return [f, self.one().value, self.zero().value] # changed

        [r0, r1] = [f, g]
        [s0, s1] = [self.one().value, self.zero().value]
        [t0, t1] = [self.zero().value, self.one().value]
        [q,r] = self.div_mod(r0, r1)
        while not(self.is_zero(r)):
            [r0, r1, s0, s1, t0, t1] = [r1, self.sub(r0, self.mult(q, r1)),
                                        s1, self.sub(s0, self.mult(q, s1)),
                                        t1, self.sub(t0, self.mult(q, t1))]
            [q,r] = self.div_mod(r0, r1)
        l = r1[-1]
        r1 = list(map(lambda x: x.div(l), r1))
        s1 = list(map(lambda x: x.div(l), s1))
        t1 = list(map(lambda x: x.div(l), t1))
        return [r1, s1, t1]

    def random_element(self, d):
        """Return a random element of degree less than d. Relies on basering having a random element function."""
        return self([self.basering.random_element() for i in range(d)])


class polynomialring_over_field_element(polynomialring_element):
    # using generic methods mostly!
    def div_mod(self, q):
        h = self.ring.div_mod(self.value, q.value)
        return [self.ring.element(h[0], self.ring), self.ring.element(h[1], self.ring)]

    def div(self, q):
        return self.div_mod(q)[0]

    def mod(self, q):
        return self.div_mod(q)[1]

    def Bezout(self, g):
        [d, a, b] = self.ring.Bezout(self.value, g.value)
        return [self.__class__(d, self.ring), self.__class__(a, self.ring), self.__class__(b, self.ring)]

    def  gcd(self, g):
         return (self.Bezout(g))[0]

    def __mod__(self, g): # overloading "%"
        return self.mod(g)

    def is_divisible_by(self, g):
        return self.ring.is_zero(self.mod(g))

    def is_irreducible(self):
    # uses Rabin irreducibility test
    # https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Rabin.27s_test_of_irreducibility
    # f must be an element of polynomialring_over_field_element, with finite coefficient field F_q
    # returns True if f is irreducible over F_q
    # returns False otherwise

        P = self.ring
        F = P.basering
        q = F.cardinality
        n = self.deg()

        def compute(e):
            r = P.one()
            b = P([0,1])
            for i in reversed(e.digits(2)):
                if i == '1':
                    r = r*b
                b = b*b
                if b.deg() > n:
                    b = b % self
            return r-P([0,1])
        
        # Step 1: The "other" checks
        one = P.one()
        h = [int(n/d) for d in [*set(pf.primefac(n))]]
        for nn in h:
            p = compute(q**nn)
            if self.gcd(p) != one:
                return False

        # Step 2: check if f divides x^q^n - x. If not: reducible.
        p = compute(q**n)
        return p.mod(self).is_zero()

        # # Step 1: check if f divides x^q^n - x. If not: reducible.
        # pp = [F.zero() for i in range(1+q**n)]
        # pp[-1] = F.one()
        # pp[1]  = -F.one()
        # p = P(pp)
        # if not(p.mod(self).is_zero()):
        #     return False

        # # Step 2: the other checks
        # one = P.one()
        # h = [int(n/d) for d in [*set(pf.primefac(n))]]
        # for nn in h:
        #     pp = [F.zero() for i in range(1+q**nn)]
        #     pp[-1] = F.one()
        #     pp[1]  = -F.one()            
        #     p = P(pp)
        #     # if not(self.gcd(p)-one).is_zero():
        #     if self.gcd(p) != one:
        #         return False

        # return True
    

class Galoisfield(field):

    def __init__(self, p, print_modulus = True): # p is a monic irreducible polynomial with prime field coefficients of degree > 0.
        # TODO: check deg(p) > 0, monic, irreducible, and so on.
        if not(isinstance(p,  polynomialring_over_field_element)):
            print("Class Galoisfield: Modulus must be a polynomial with coefficients in a finite prime field")
            raise
        if not(isinstance(p.ring.basering, primefield)):
            print("Class Galoisfield: Modulus polynomial must have coefficients in a finite prime field")
            raise
        self.print_modulus = print_modulus
        self.modulus = p
        self.p_ring = p.ring                          # the underlying polynomial ring
        self.c_ring = self.p_ring.basering            # coefficients of polynomials
        self.char = self.c_ring.modulus
        self.deg = p.deg()
        self.cardinality = self.char ** self.deg
        self.iterator = globals()['Galoisfield_iterator']
        super().__init__()

    def __str__(self):
        return 'GF(' + str(self.char) + '^' + str(self.deg) + ')'

    def __repr__(self):
        return self.__class__.__name__ + '(' + self.modulus.__repr__() + ')'

    def normalise(self, v):  # v is an object of type self.p_ring, that is, a polynomial with coefficients in finite prime field
        return v.mod(self.modulus)

    def zero(self):
        return self.element(self.p_ring.zero(), self)

    def one(self):
        return self.element(self.p_ring.one(), self)

    def add(self, a, b):
        return a.add(b)

    def sub(self, a, b):
        return a.sub(b)

    def mult(self, a, b):
        return self.normalise(a.mult(b))

    def div(self, a, b): # compute a/b
        if b.is_zero():
            print("Class Galoisfield: Cannot divide by 0")
            raise
        return a.mult(b.Bezout(self.modulus)[1])

    def is_zero(self, v):
        return v.is_zero()

    def __iter__(self):
        return self.iterator(self)

    def random_element(self):
        """Return a random element."""
        return self([self.c_ring.random_element() for i in range(self.deg)])


class Galoisfield_element(field_element):
    # using generic methods mostly!
    def __str__(self):
        if self.ring.print_modulus:
            return str(self.value) + ' mod ' + str(self.ring.modulus)
        else:
            return str(self.value)

    def __init__(self, v, G):
        if not(isinstance(G, Galoisfield)):
            print("Class Galoisfield_element: Specify a Galois field")
            raise
        w = v
        if not(isinstance(v, polynomialring_over_field_element)):
            w = polynomialring_over_field_element(v, G.p_ring)
        super().__init__(w, G)


class Galoisfield_iterator:
    def __init__(self, G):
        self.current = 0
        self.max     = G.char**G.deg - 1
        self.field   = G

    def __next__(self):
        n = self.current
        if n > self.max:
            raise StopIteration
        nums = []
        for i in range(self.field.deg):
            n, r = divmod(n, self.field.char)
            nums.append(str(r))
        self.current += 1
        return self.field.element(self.field.p_ring.element(nums, self.field.p_ring), self.field)
# Rings and fields:1 ends here
