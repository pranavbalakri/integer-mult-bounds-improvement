"""Checks of the new guard, distinct from the inherited network audit."""
import itertools
import unittest
from fractions import Fraction as Q
from semantic_guard import certify


def mul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def conj(a): return a[0], -a[1]
def add(a, b): return a[0]+b[0], a[1]+b[1]


def phase_numerator(e, phase):
    """Integer numerator of H diag(i^phase) H, denominator 2^e."""
    roots = [(1,0),(0,1),(-1,0),(0,-1)]
    n = 1 << e
    out = []
    for x in range(n):
        row = []
        for y in range(n):
            a = (0,0)
            for z in range(n):
                sign = -1 if bin((x^y)&z).count('1') % 2 else 1
                r = roots[phase[z] % 4]
                a = add(a, (sign*r[0],sign*r[1]))
            row.append(a)
        out.append(row)
    return out


class SemanticGuard(unittest.TestCase):
    def test_every_two_bit_phase(self):
        # Exhaust all 4^4 phase functions, not just quadratic/example frames.
        # Numerators are Gaussian integers; row orthogonality certifies the
        # exact unitarity used in the infinity-norm Cauchy--Schwarz bound.
        e = 2; n = 1 << e
        for p in itertools.product(range(4), repeat=n):
            a = phase_numerator(e, p)
            for i in range(n):
                for j in range(n):
                    total = (0,0)
                    for k in range(n): total = add(total, mul(a[i][k],conj(a[j][k])))
                    self.assertEqual(total, (n*n if i == j else 0,0))

    def test_semantic_return_cancels_denominators(self):
        # C^4=I: a computation-tree depth count grows, but the returned value
        # has exactly its original denominator. This is why sibling guards
        # may not be summed as if every internal halving survived its return.
        C = [[(Q(1,2),Q(1,2)),(Q(1,2),Q(-1,2))],
             [(Q(1,2),Q(-1,2)),(Q(1,2),Q(1,2))]]
        v = [(Q(1,8),Q(3,16)),(Q(-5,32),Q(1,4))]; original = list(v)
        for _ in range(4):
            v = [add(mul(row[0],v[0]),mul(row[1],v[1])) for row in C]
        self.assertEqual(v, original)

    def test_h16_scalar_and_guard_certificate(self):
        c = certify(16)
        self.assertEqual(c['forward_G'], '135386')
        self.assertEqual(c['inverse_G'], '256057')
        self.assertEqual(c['G'], '69333066004')
        self.assertEqual(c['B'], 2)
        self.assertEqual(c['recursive_guard_coefficient'], 36)
        self.assertEqual(c['C1'], 1)
        self.assertEqual(c['precision_exponent'], 3)

    def test_cubic_numeric_margins(self):
        # Stronger uniform check d<=b and L^eps<=b. This isolates the numeric
        # inequalities from the inherited eventual prime/layout hypotheses.
        for b in range(1, 101):
            for d in range(1, b+1):
                p = 768*b**3
                self.assertLessEqual(128*(d+1), p)
                alpha2 = (16*b**3)//d
                self.assertGreaterEqual(Q(alpha2,8*d*b), Q(15,8))
                self.assertLessEqual(2*d*alpha2, 32*b**3)

    def test_integer_remainders(self):
        m=256; h=16; C=36
        for e in range(m, 5*m+1):
            r=e % m; b=e-r
            if b:
                self.assertLessEqual(2*b+64+C*(m-h)*(b//m),C*b)
                self.assertLessEqual(r+C*b,C*e)
            self.assertLessEqual(8*r,C*e)
        # A smaller unsupported slope is rejected at the smallest internal node.
        self.assertGreater(2*m+64+35*(m-h),35*m)


if __name__ == '__main__': unittest.main()
