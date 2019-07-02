#!/usr/bin/python

"""Pairings over a 256-bit BN curve
(C) 2017 Jack Lloyd <jack@randombit.net>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import hashlib
import sys
import os

v = 1868033
#v = 0b111001000000100000001
u = pow(v, 3)

#p = 36*pow(u,4) + 36*pow(u,3) + 24*pow(u,2) + 6*u + 1
#order = 36*pow(u,4) + 36*pow(u,3) + 18*pow(u,2) + 6*u + 1

p = (((u + 1)*6*u + 4)*u + 1)*6*u + 1
order = p - 6*u*u

def is_py3():
    return (sys.version_info[0] == 3)

def is_integer_type(x):
    if is_py3():
        return type(x) in [int]
    else:
        return type(x) in [int,long]

def inverse_mod(a, n):
    t = 0
    t2 = 1
    r = n
    r2 = a

    while r2 != 0:
        q = r // r2
        (t, t2) = (t2, t - q * t2)
        (r, r2) = (r2, r - q * r2)

    if r > 1:
        return 0
    if t < 0:
        t += n
    return t

def sqrt_mod_p(a):
    assert p % 4 == 3
    return pow(a, (p+1)//4, p)

def inv_mod_p(a):
    # Fermat
    return pow(a, p-2, p)

def legendre(a):
    x = pow(a, (p-1)//2, p)
    if x == 0 or x == 1:
        return x
    if x == p-1:
        return -1
    assert False

# Takes larger random range to minimize ratio of unusable numbers
rand_elem_bytes = (order.bit_length() + 7) // 8 + 1
rand_elem_base = 2
rand_elem_range = order - rand_elem_base
rand_elem_barrier = (1 << (8 * rand_elem_bytes)) - rand_elem_range

def rand_elem():
    """ Debiased random element generator """
    while True:
        rand_bytes = os.urandom(rand_elem_bytes)
        rand_num = int.from_bytes(rand_bytes, sys.byteorder)
        res = rand_num % rand_elem_range
        if (rand_num - res) <= rand_elem_barrier:
            return res + rand_elem_base

# Montgomery params
R = pow(2,256)
R1 = R % p
R2 = (R*R) % p
R3 = (R1*R2) % p
N = R - inverse_mod(p, R)

def to_naf(x):
    z = []
    while x > 0:
        if x % 2 == 0:
            z.append(0)
        else:
            zi = 2 - (x % 4)
            x -= zi
            z.append(zi)
        x = x // 2
    return z

# 6u+2 in NAF
naf_6up2 = list(reversed(to_naf(6*u+2)))[1:]

def bits_of(k):
    return [int(c) for c in "{0:b}".format(k)]

class gfp_1(object):
    def _redc(self, T):
        #assert T < (R*p-1)
        m = ((T & (R-1)) * N) & (R-1)
        t = (T + m*p) >> 256
        if t >= p:
            t -= p
        return t

    def __init__(self, v, redc_needed = True):
        #assert v >= 0
        if redc_needed:
            #assert v < R
            self.v = self._redc(v * R2)
            #assert self.value() == v % p
        else:
            #assert v < p
            self.v = v

    def __eq__(self, other):
        return self.v == other.v

    def __ne__(self, other):
        return self.v != other.v

    def __str__(self):
        return "%d" % (self.value())

    def __add__(self, other):
        x = self.v + other.v
        if (x >> 255) > 0 and x >= p:
            x -= p
        #assert self._redc(x) == (self.value() + other.value()) % p
        return gfp_1(x, False)

    def __sub__(self,other):
        x = self.v - other.v
        if x < 0:
            x += p
        #assert self._redc(x) == (self.value() - other.value()) % p
        return gfp_1(x, False)

    def __mul__(self,other):
        return gfp_1(self._redc(self.v * other.v), False)

    def square(self):
        return self * self
        return gfp_1(self.v * self.v)

    def double(self):
        return self + self

    def triple(self):
        return self + self + self

    def is_one(self):
        return self.v == R1

    def is_zero(self):
        return self.v == 0

    def value(self):
        return self._redc(self.v)

    def inverse(self):
        # Fermat
        x = gfp_1(self._redc(R3 * inv_mod_p(self.v)), False)
        #assert (self * x).value() == 1
        return x

    def additive_inverse(self):
        x = gfp_1(p, True) - self
        #assert (self + x).value() == 0
        return x

    def to_bytes(self):
        p_bytes = (p.bit_length() // 8) + (1 if (p.bit_length() % 8) > 0 else 0)
        return self.value().to_bytes(p_bytes, 'big')

def point_on_curve(point, b):
    point.force_affine()
    yy = point.y.square()
    xxx = point.x.square() * point.x
    yy -= xxx
    yy -= b
    return yy.is_zero()

def point_add(a, b):
    if a.is_infinite():
        return b
    if b.is_infinite():
        return a

    """
    http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
      Z1Z1 = a.z^2
      Z2Z2 = b.z^2
      U1 = a.x*Z2Z2
      U2 = b.x*Z1Z1
      S1 = a.y*b.z*Z2Z2
      S2 = b.y*a.z*Z1Z1
      H = U2-U1
      I = (2*H)^2
      J = H*I
      r = 2*(S2-S1)
      V = U1*I
      X3 = r^2-J-2*V
      Y3 = r*(V-X3)-2*S1*J
      Z3 = ((a.z+b.z)^2-Z1Z1-Z2Z2)*H
        """

    z1z1 = a.z.square()
    z2z2 = b.z.square()
    u1 = (z2z2 * a.x)
    u2 = (z1z1 * b.x)
    h = u2 - u1

    s1 = (a.y * b.z * z2z2)
    s2 = (b.y * a.z * z1z1)
    r = s2 - s1

    if h.is_zero() and r.is_zero():
        return a.double()

    r = r.double()
    i = h.square()
    i = i.double().double()
    j = (h * i)

    V = (u1 * i)

    c_x = (r.square() - j - V.double())
    c_y = (r * (V - c_x) - s1*j.double())

    c_z = a.z + b.z
    c_z = c_z.square()
    c_z -= z1z1
    c_z -= z2z2
    c_z *= h

    return a.__class__(c_x, c_y, c_z)

def point_double(a):
    # http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
    """
   
    compute A = X1^2
        compute B = Y1^2
        compute C = B^2
        compute D = 2 ((X1 + B)^2 - A - C)
        compute E = 3 A
        compute F = E^2
        compute X3 = F - 2 D
        compute Y3 = E (D - X3) - 8 C
        compute Z3 = 2 Y1 Z1
    """
    A = a.x.square()
    B = a.y.square()
    C = B.square()

    t = a.x + B
    t = t.square()

    D = (t - A - C)
    D = D.double()

    E = A.double() + A
    F = E.square()

    C8 = C.double().double().double()

    c_x = (F - D.double())
    c_y = (E * (D - c_x) - C8)
    c_z = (a.y * a.z).double()

    return a.__class__(c_x, c_y, c_z)

def point_force_affine(point):
    if point.z.is_one():
        return

    zinv = point.z.inverse()
    zinv2 = (zinv * zinv)
    zinv3 = (zinv2 * zinv)

    point.x = point.x * zinv2
    point.y = point.y * zinv3
    point.z = point.one_element()

def point_scalar_mul(pt, k):
    assert is_integer_type(k)

    R = [pt.__class__(pt.zero_element(), pt.zero_element(), pt.zero_element()),
         pt]

    for kb in bits_of(k):
        R[kb^1] = R[kb].add(R[kb^1])
        R[kb] = R[kb].double()
    return R[0]

curve_B = gfp_1(3)

class curve_point(object):
    def __init__(self, x, y, z = gfp_1(1)):
        assert type(x) in [gfp_1]
        assert type(y) in [gfp_1]
        assert type(z) in [gfp_1]

        self.x = x
        self.y = y
        self.z = z

    def zero_element(self):
        return gfp_1(0)

    def one_element(self):
        return gfp_1(1)

    def __repr__(self):
        self.force_affine()
        return "(%d, %d)" % (self.x.value(), self.y.value())

    def is_on_curve(self):
        return point_on_curve(self, curve_B)

    def is_infinite(self):
        return self.z.is_zero()

    def add(a, b):
        return point_add(a, b)

    def double(a):
        return point_double(a)

    def force_affine(self):
        point_force_affine(self)

    def scalar_mul(self, k):
        return point_scalar_mul(self, k)

# Any point (1,y) where y is a square root of b+1 is a generator
curve_G = curve_point(gfp_1(1), gfp_1(p-2))

assert curve_G.is_on_curve()

class gfp_2(object):
    def __init__(self, x, y):
        """
        Represented as i*x + y
        """
        if type(x) == gfp_1:
            self.x = x
            self.y = y
        else:
            # Assumed to be integers
            self.x = gfp_1(x)
            self.y = gfp_1(y)

    def __repr__(self):
        return "(%d,%d)" % (self.x.value(), self.y.value())

    def __eq__(self,other):
        return self.x == other.x and self.y == other.y

    def __ne__(self,other):
        return self.x != other.x or self.y != other.y

    def is_zero(self):
        return self.x.is_zero() and self.y.is_zero()

    def is_one(self):
        return self.x.is_zero() and self.y.is_one()

    def conjugate_of(self):
        """
        For gamma = A + iB \in gfp2
        gamma^p = A - iB
        """
        return gfp_2(self.x.additive_inverse(), self.y)

    def negative_of(self):
        return gfp_2(self.x.additive_inverse(), self.y.additive_inverse())

    def add(self, other):
        return gfp_2((self.x + other.x), (self.y + other.y))

    def sub(self, other):
        return gfp_2((self.x - other.x), (self.y - other.y))

    def double(self):
        return gfp_2((self.x.double()), (self.y.double()))

    def mul(a, b):
        assert type(a) == gfp_2 and type(b) == gfp_2
        # Karatsuba
        vy = (a.y * b.y)
        vx = (a.x * b.x)
        c0 = (vy - vx)
        c1 = ((a.x + a.y)*(b.x + b.y) - vy - vx)

        return gfp_2(c1,c0)

    def __mul__(a,b):
        return a.mul(b)

    def __sub__(a,b):
        return a.sub(b)

    def __add__(a,b):
        return a.add(b)

    def mul_scalar(self, k):
        return gfp_2((self.x * k), (self.y * k))

    # Multiply by i+3
    def mul_xi(a):
        # (xi + y)(3 + i) = 3xi + 3y - x + yi = (3x + y)i + (3y - x)
        tx = (a.x.triple()) + a.y
        ty = (a.y.triple()) - a.x
        return gfp_2(tx, ty)

    def square(a):
        assert type(a.x) == gfp_1
        assert type(a.y) == gfp_1
        # Complex squaring
        t1 = a.y - a.x
        t2 = a.y + a.x
        ty = (t1 * t2)
        #ty = a.y*a.y - a.x*a.x
        tx = (a.x * a.y)
        tx = tx.double()
        return gfp_2(tx, ty)

    def inverse(a):
        # Algorithm 8 from http://eprint.iacr.org/2010/354.pdf
        t = a.x.square() + a.y.square()

        inv = t.inverse()

        c_x = (a.x.additive_inverse() * inv)
        c_y = (a.y * inv)

        return gfp_2(c_x, c_y)

    def exp(p, k):
        assert is_integer_type(k)
        assert type(p) == gfp_2
        R = [gfp_2(gfp_1(0),gfp_1(1)),p]
        for kb in bits_of(k):
            R[kb^1] = R[kb].mul(R[kb^1])
            assert type(R[kb]) == gfp_2
            R[kb] = R[kb].square()
        return R[0]

gfp_2_zero = gfp_2(0, 0)
gfp_2_one = gfp_2(0, 1)

xi = gfp_2(1,3) # i + 3

xi1 = [
    xi.exp(1*(p-1)//6),
    xi.exp(2*(p-1)//6),
    xi.exp(3*(p-1)//6),
    xi.exp(4*(p-1)//6),
    xi.exp(5*(p-1)//6)
]

xi2 = [(x * x.conjugate_of()) for x in xi1]

# twist of the curve over GF(p^2)

twist_B = xi.inverse().mul(gfp_2(0, curve_B.value()))

class curve_twist(object):
    def __init__(self, x, y, z):
        assert type(x) == gfp_2 and type(y) == gfp_2 and type(z) == gfp_2
        self.x = x
        self.y = y
        self.z = z

    def one_element(self):
        return gfp_2_one

    def zero_element(self):
        return gfp_2_zero

    def __repr__(self):
        self.force_affine()
        return "(%s, %s)" % (self.x, self.y)

    def is_on_curve(self):
        return point_on_curve(self, twist_B)

    def is_infinite(self):
        return self.z.is_zero()

    # Add two points on the twist
    def add(a, b):
        return point_add(a,b)

    def double(a):
        return point_double(a)

    def scalar_mul(self, k):
        return point_scalar_mul(self, k)

    def force_affine(self):
        point_force_affine(self)

    def negate(self):
        self.y = self.y.negative_of()

# TODO derive this
twist_G = curve_twist(
    gfp_2(21167961636542580255011770066570541300993051739349375019639421053990175267184,
          64746500191241794695844075326670126197795977525365406531717464316923369116492),
    gfp_2(20666913350058776956210519119118544732556678129809273996262322366050359951122,
          17778617556404439934652658462602675281523610326338642107814333856843981424549),
    gfp_2(0,1))

assert twist_G.is_on_curve()

# cubic extension of gfp_2
class gfp_6(object):
    def __init__(self, x, y, z):
        assert type(x) == gfp_2 and type(y) == gfp_2 and type(z) == gfp_2
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self,other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __repr__(self):
        return "(%s,%s,%s)" % (self.x, self.y, self.z)

    def is_zero(self):
        return self.x.is_zero() and self.y.is_zero() and self.z.is_zero()

    def is_one():
        return self.x.is_zero() and self.y.is_zero() and self.z.is_one()

    def negative_of(self):
        return gfp_6(self.x.negative_of(), self.y.negative_of(), self.z.negative_of())

    def add(a, b):
        return gfp_6(a.x.add(b.x),
                     a.y.add(b.y),
                     a.z.add(b.z))

    def sub(a, b):
        return gfp_6(a.x.sub(b.x),
                     a.y.sub(b.y),
                     a.z.sub(b.z))

    def double(self):
        return gfp_6(self.x.double(),
                     self.y.double(),
                     self.z.double())

    def mul(a, b):
        # Algorithm 13 from http://eprint.iacr.org/2010/354.pdf
        # plus some short-circuits

        if a.x.is_zero():
            if a.y.is_zero():
                return b.mul_scalar(a.z)

            t0 = (b.z * a.z)
            t1 = (b.y * a.y)

            tz = (b.x + b.y) * (a.y)
            tz -= t1
            tz = tz.mul_xi()
            tz += t0

            ty = (b.y + b.z) * (a.y + a.z)
            ty -= t0
            ty -= t1

            tx = (b.x) * (a.z)
            tx += t1

            return gfp_6(tx, ty, tz)

        if b.x.is_zero():
            if b.y.is_zero():
                return a.mul_scalar(b.z)

            t0 = (a.z * b.z)
            t1 = (a.y * b.y)

            tz = (a.x + a.y) * (b.y)
            tz -= t1
            tz = tz.mul_xi()
            tz += t0

            ty = (a.y + a.z) * (b.y + b.z)
            ty -= t0
            ty -= t1

            tx = (a.x) * (b.z)
            tx += t1

            return gfp_6(tx, ty, tz)

        t0 = (a.z * b.z)
        t1 = (a.y * b.y)
        t2 = (a.x * b.x)

        tz = (a.x + a.y) * (b.x + b.y)
        tz -= t1
        tz -= t2
        tz = tz.mul_xi()
        tz += t0

        ty = (a.y + a.z) * (b.y + b.z)
        ty -= t0
        ty -= t1
        ty += t2.mul_xi()

        tx = (a.x + a.z) * (b.x + b.z)
        tx -= t0
        tx += t1
        tx -= t2

        return gfp_6(tx, ty, tz)

    def __mul__(a,b):
        return a.mul(b)

    def __add__(a,b):
        return a.add(b)
    def __sub__(a,b):
        return a.sub(b)

    def mul_scalar(self, k):
        assert type(k) == gfp_2

        return gfp_6(self.x.mul(k),
                     self.y.mul(k),
                     self.z.mul(k))

    def mul_tau(a):
        tx = a.y
        ty = a.z
        tz = a.x.mul_xi()
        return gfp_6(tx, ty, tz)

    def square(a):
        # Algorithm 16 from http://eprint.iacr.org/2010/354.pdf
        ay2 = a.y.double()
        c4 = (a.z * ay2)
        c5 = a.x.square()
        c1 = c5.mul_xi() + c4
        c2 = c4 - c5
        c3 = a.z.square()
        c4 = a.x + a.z - a.y
        c5 = (ay2 * a.x)
        c4 = c4.square()
        c0 = c5.mul_xi() + c3
        c2 = c2 + c4 + c5 - c3
        n = gfp_6(c2, c1, c0)
        return n

    def inverse(a):
        # Algorithm 17
        XX = a.x.square()
        YY = a.y.square()
        ZZ = a.z.square()

        XY = (a.x * a.y)
        XZ = (a.x * a.z)
        YZ = (a.y * a.z)

        A = ZZ - XY.mul_xi()
        B = XX.mul_xi() - YZ
        # There is an error in the paper for this line
        C = YY - XZ

        F = (C * a.y).mul_xi()
        F += (A * a.z)
        F += (B * a.x).mul_xi()

        F = F.inverse()

        c_x = C * F
        c_y = B * F
        c_z = A * F
        return gfp_6(c_x, c_y, c_z)

gfp_6_zero = gfp_6(gfp_2_zero, gfp_2_zero, gfp_2_zero)
gfp_6_one  = gfp_6(gfp_2_zero, gfp_2_zero, gfp_2_one)

class gfp_12(object):
    def __init__(self, x, y = None):
        assert type(x) == gfp_6
        assert type(y) == gfp_6
        self.x = x
        self.y = y

    def __eq__(self,other):
        return self.x == other.x and self.y == other.y

    def __repr__(self):
        return "(%s,%s)" % (self.x, self.y)

    def is_zero(self):
        return self.x.is_zero() and self.y.is_zero()

    def is_one(self):
        return self.x.is_zero() and self.y.is_one()

    def conjugate_of(self):
        return gfp_12(self.x.negative_of(), self.y)

    def negative_of(self):
        return gfp_12(self.x.negative_of(), self.y.negative_of())

    def frobenius(self):
        e1_x = self.x.x.conjugate_of().mul(xi1[4])
        e1_y = self.x.y.conjugate_of().mul(xi1[2])
        e1_z = self.x.z.conjugate_of().mul(xi1[0])

        e2_x = self.y.x.conjugate_of().mul(xi1[3])
        e2_y = self.y.y.conjugate_of().mul(xi1[1])
        e2_z = self.y.z.conjugate_of()

        return gfp_12(gfp_6(e1_x,e1_y,e1_z), gfp_6(e2_x,e2_y,e2_z))

    def frobenius_p2(self):
        e1_x = self.x.x.mul(xi2[4])
        e1_y = self.x.y.mul(xi2[2])
        e1_z = self.x.z.mul(xi2[0])

        e2_x = self.y.x.mul(xi2[3])
        e2_y = self.y.y.mul(xi2[1])
        e2_z = self.y.z

        return gfp_12(gfp_6(e1_x,e1_y,e1_z), gfp_6(e2_x,e2_y,e2_z))

    def sub(a, b):
        return gfp_12(a.x - b.x, a.y - b.y)

    def mul(a, b):
        # TODO Karatsuba (algo 20)
        AXBX = a.x * b.x
        AXBY = a.x * b.y
        AYBX = a.y * b.x
        AYBY = a.y * b.y
        return gfp_12(AXBY + AYBX, AYBY + AXBX.mul_tau())

    def mul_scalar(self, k):
        assert type(k) == gfp_6
        return gfp_12(self.x.mul(k), self.y.mul(k))

    def exp(self, k):
        assert is_integer_type(k)

        R = [gfp_12(gfp_6_zero, gfp_6_one), self]

        for kb in bits_of(k):
            R[kb^1] = R[kb].mul(R[kb^1])
            R[kb] = R[kb].square()

        return R[0]

    def square(a):
        v0 = a.x * a.y
        t = a.x.mul_tau()
        t += a.y
        ty = a.x + a.y
        ty *= t
        ty -= v0
        t = v0.mul_tau()
        ty -= t

        c_x = v0.double()
        c_y = ty

        return gfp_12(c_x, c_y)

    def inverse(a):
        e = gfp_12(a.x.negative_of(), a.y)

        t1 = a.x.square()
        t2 = a.y.square()
        t1 = t1.mul_tau()
        t1 = t2 - t1
        t2 = t1.inverse()

        e = e.mul_scalar(t2)
        return e


def line_func_add(r, p, q, r2):
    assert type(r) == curve_twist
    assert type(p) == curve_twist
    assert type(q) == curve_point
    assert type(r2) == gfp_2

    r_t = r.z.square()
    B = p.x * r_t
    D = p.y + r.z
    D = D.square()
    D -= r2
    D -= r_t
    D *= r_t

    H = B - r.x
    I = H.square()

    E = I.double().double()

    J = H * E
    L1 = D - r.y
    L1 -= r.y

    V = r.x * E

    r_x = L1.square()
    r_x -= J
    r_x -= V.double()

    r_z = r.z + H
    r_z = r_z.square()
    r_z -= r_t
    r_z -= I

    t = V - r_x
    t *= L1
    t2 = r.y * J
    t2 = t2.double()
    r_y = t - t2

    r_out = curve_twist(r_x, r_y, r_z)

    t = p.y + r_z
    t = t.square()
    t = t - r2
    t = t - (r_z.square())

    t2 = L1 * p.x
    t2 = t2.double()
    a = t2 - t

    c = r_z.mul_scalar(q.y).double()

    b = L1.negative_of()
    b = b.mul_scalar(q.x).double()

    return (a, b, c, r_out)

def line_func_double(r, q):
    assert type(r) == curve_twist
    assert type(q) == curve_point

    # cache this?
    r_t = r.z.square()

    A = r.x.square()
    B = r.y.square()
    C = B.square()

    D = r.x + B
    D = D.square()
    D -= A
    D -= C
    D = D.double()

    E = A.double() + A
    F = E.square()

    C8 = C.double().double().double() # C*8

    r_x = F - D.double()
    r_y = E * (D - r_x) - C8

    # (y+z)*(y+z) - (y*y) - (z*z) = 2*y*z
    r_z = (r.y + r.z).square() - B - r_t

    assert r_z == r.y*r.z.double()

    r_out = curve_twist(r_x, r_y,r_z)
    #assert r_out.is_on_curve()

    a = r.x + E
    a = a.square()
    a -= (A + F + B.double().double())

    t = E * r_t
    t = t.double()
    b = t.negative_of()
    b = b.mul_scalar(q.x)

    c = r_z * r_t
    c = c.double().mul_scalar(q.y)

    return (a,b,c,r_out)

def mul_line(r, a, b, c):
    assert type(r) == gfp_12
    assert type(a) == gfp_2
    assert type(b) == gfp_2
    assert type(c) == gfp_2

    # See function fp12e_mul_line in dclxvi

    t1 = gfp_6(gfp_2_zero, a, b)
    t2 = gfp_6(gfp_2_zero, a, b + c)

    t1 = t1 * r.x
    t3 = r.y.mul_scalar(c)
    r.x += r.y
    r.y = t3
    r.x *= t2
    r.x -= t1
    r.x -= r.y
    r.y += t1.mul_tau()

def miller(q, p):

    import copy

    assert type(q) == curve_twist
    assert type(p) == curve_point

    Q = copy.deepcopy(q)
    Q.force_affine()

    P = copy.deepcopy(p)
    P.force_affine()

    mQ = copy.deepcopy(Q)
    mQ.negate()

    f = gfp_12(gfp_6_zero, gfp_6_one)
    T = Q

    Qp = Q.y.square()

    for naf_i in naf_6up2:
        # Skip on first iteration?
        f = f.square()

        a, b, c, T = line_func_double(T, P)
        mul_line(f, a, b, c)

        if naf_i == 1:
            a, b, c, T = line_func_add(T, Q, P, Qp)
            mul_line(f, a, b, c)
        elif naf_i == -1:
            a, b, c, T = line_func_add(T, mQ, P, Qp)
            mul_line(f, a, b, c)

    # Q1 = pi(Q)
    Q1 = curve_twist(
        Q.x.conjugate_of().mul(xi1[1]),
        Q.y.conjugate_of().mul(xi1[2]),
        gfp_2_one)

    # Q2 = pi2(Q)
    Q2 = curve_twist(
        Q.x.mul_scalar(xi2[1].y),
        Q.y,
        gfp_2_one)

    Qp = Q1.y.square()
    a, b, c, T = line_func_add(T, Q1, P, Qp)
    mul_line(f, a, b, c)

    Qp = Q2.y.square()
    a, b, c, T = line_func_add(T, Q2, P, Qp)
    mul_line(f, a, b, c)

    return f

def final_exp(inp):
    assert type(inp) == gfp_12

    # Algorithm 31 from https://eprint.iacr.org/2010/354.pdf

    t1 = inp.conjugate_of()
    inv = inp.inverse()

    t1 = t1.mul(inv)
    # Now t1 = inp^(p**6-1)

    t2 = t1.frobenius_p2()
    t1 = t1.mul(t2)

    fp1 = t1.frobenius()
    fp2 = t1.frobenius_p2()
    fp3 = fp2.frobenius()

    fu1 = t1.exp(u)
    fu2 = fu1.exp(u)
    fu3 = fu2.exp(u)

    y3 = fu1.frobenius()
    fu2p = fu2.frobenius()
    fu3p = fu3.frobenius()
    y2 = fu2.frobenius_p2()

    y0 = fp1.mul(fp2)
    y0 = y0.mul(fp3)

    y1 = t1.conjugate_of()
    y5 = fu2.conjugate_of()
    y3 = y3.conjugate_of()
    y4 = fu1.mul(fu2p)
    y4 = y4.conjugate_of()

    y6 = fu3.mul(fu3p)
    y6 = y6.conjugate_of()

    t0 = y6.square()
    t0 = t0.mul(y4)
    t0 = t0.mul(y5)

    t1 = y3.mul(y5)
    t1 = t1.mul(t0)
    t0 = t0.mul(y2)
    t1 = t1.square()
    t1 = t1.mul(t0)
    t1 = t1.square()
    t0 = t1.mul(y1)
    t1 = t1.mul(y0)
    t0 = t0.square()
    t0 = t0.mul(t1)

    return t0

def optimal_ate(a, b):
    assert type(a) == curve_twist
    assert type(b) == curve_point

    e = miller(a, b)
    ret = final_exp(e)

    if a.is_infinite() or b.is_infinite():
        return gfp_12(gfp_6_zero, gfp_6_one)

    return ret

def g1_scalar_base_mult(k):
    return curve_G.scalar_mul(k)

def g1_random():
    k = rand_elem()
    return k, g1_scalar_base_mult(k)

def g1_add(a, b):
    assert type(a) == curve_point
    assert type(b) == curve_point

    return a.add(b)

def g1_marshall(a):
    a.force_affine()
    return (a.x,a.y)

def g1_unmarshall(x,y):
    return curve_point(x, y)

def g1_hash_to_point(msg):
    # From "Indifferentiable Hashing to Barreto-Naehrig Curves"
    # https://www.di.ens.fr/~fouque/pub/latincrypt12.pdf

    # constants
    sqrt_neg_3 = sqrt_mod_p(p-3)
    inv_2 = inv_mod_p(2)
    b = curve_B.value()

    # compute t in F_q
    sha = hashlib.sha512()
    sha.update(msg)
    t = int(sha.hexdigest(), 16) % p

    if t == 0:
        # TODO handle this case as described in paper
        assert False

    t2 = (t*t) % p

    chi_t = legendre(t)

    w = sqrt_neg_3 * t * inv_mod_p(1 + b + t2)

    def g(x):
        return (x*x*x + b) % p

    x1 = ((sqrt_neg_3 - 1) * inv_2 - t*w) % p
    g_x1 = g(x1)
    if legendre(g_x1) == 1:
        x1_sqrt = sqrt_mod_p(g_x1)
        return curve_point(gfp_1(x1),
                           gfp_1(chi_t * x1_sqrt))

    x2 = (-1 - x1) % p
    g_x2 = g(x2)

    if legendre(g_x2) == 1:
        x2_sqrt = sqrt_mod_p(g_x2)
        return curve_point(gfp_1(x2),
                           gfp_1(chi_t * x2_sqrt))

    x3 = 1 + inv_mod_p(w*w)
    g_x3 = g(x3)

    assert legendre(g_x3) == 1
    x3_sqrt = sqrt_mod_p(g_x3)
    return curve_point(gfp_1(x3),
                       gfp_1(chi_t * x3_sqrt))

def g1_compress(g1):
    g1.force_affine()
    x = g1.x.value()
    y = g1.y.value()

    return (x, y & 1)

def g1_uncompress(g):
    x = g[0]
    y_sign = g[1]

    assert y_sign == 0 or y_sign == 1
    assert x >= 0 and x < p

    xxx = (x*x*x + curve_B.value()) % p

    y = sqrt_mod_p(xxx)

    if y_sign != y & 1:
        y = p - y

    return curve_point(gfp_1(x), gfp_1(y))

def g2_scalar_base_mult(k):
    return twist_G.scalar_mul(k)

def g2_random():
    k = rand_elem()
    return k, g2_scalar_base_mult(k)

def g2_add(a, b):
    return a.add(b)

def g2_marshall(a):
    return (a.x.x, a.x.y, a.y.x, a.y.y)

def g2_unmarshall(w,x,y,z):
    if w == x == y == z == 0:
        # This is the point at infinity.
        return curve_twist(gfp_2(0, 0), gfp_2(0, 1), gfp_2(0, 0))
    else:
        return curve_twist(gfp_2(w, x), gfp_2(y, z), gfp_2(0, 1))

def gt_scalar_mult(x, k):
    return x.exp(k)

def gt_add(a, b):
    return a.mul(b)

def gt_marshall(gt):
    return (gt.x.x.x,
            gt.x.x.y,
            gt.x.y.x,
            gt.x.y.y,
            gt.x.z.x,
            gt.x.z.y,
            gt.y.x.x,
            gt.y.x.y,
            gt.y.y.x,
            gt.y.y.y,
            gt.y.z.x,
            gt.y.z.y)

def gt_unmarshall(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11):
    return gfp_12(
        gfp_6(gfp_2(p0,p1), gfp_2(p2,p3), gfp_2(p4,p5)),
        gfp_6(gfp_2(p6,p7), gfp_2(p8,p9), gfp_2(p10,p11)))

def gt_hash(gt):
    sha = hashlib.sha512()

    for parts in gt_marshall(gt):
        parts = parts.to_bytes()
        sha.update(parts)

    return sha.digest()
