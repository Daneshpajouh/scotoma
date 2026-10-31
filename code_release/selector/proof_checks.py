#!/usr/bin/env python3
"""Exact polynomial checks for the middle-four proof, using Fraction only."""
import json
from fractions import Fraction as F


class P:
    """Four-variable polynomial over Q; variables are m,N,t,r."""
    def __init__(self, terms):
        if isinstance(terms, P):
            self.terms = dict(terms.terms)
        elif isinstance(terms, dict):
            self.terms = {k: F(v) for k,v in terms.items() if v}
        else:
            self.terms = {(0,0,0,0): F(terms)} if terms else {}

    def __add__(self, other):
        out = dict(self.terms)
        for k,v in P(other).terms.items():
            out[k] = out.get(k,F(0)) + v
        return P(out)
    __radd__ = __add__

    def __neg__(self):
        return P({k:-v for k,v in self.terms.items()})

    def __sub__(self, other):
        return self + (-P(other))

    def __rsub__(self, other):
        return P(other) + (-self)

    def __mul__(self, other):
        out = {}
        for k,v in self.terms.items():
            for q,w in P(other).terms.items():
                key = tuple(a+b for a,b in zip(k,q))
                out[key] = out.get(key,F(0)) + v*w
        return P(out)
    __rmul__ = __mul__

    def __pow__(self, n):
        out = P(1)
        for _ in range(n):
            out = out*self
        return out

    def sub(self, i, replacement):
        out = P(0)
        for k,v in self.terms.items():
            q = list(k)
            q[i] = 0
            out += P({tuple(q):v}) * replacement**k[i]
        return out

    def __eq__(self, other):
        return self.terms == P(other).terms


def var(i):
    k = [0]*4
    k[i] = 1
    return P({tuple(k):F(1)})


m,N,t,r = [var(i) for i in range(4)]
A = 9*N*N + 12*N*m - 87*N - 5*m*m - 23*m + 174
# Numerator obtained directly from the binomial ratios in the proof.
M = N-2
derived_A = (9*(M-m)*(M-m+1) + 30*(m-2)*(M-m+1)
             + 16*(m-2)*(m-3))
assert A == derived_A
H = (m-1)*A - 4*(2*N-m-3)*(N-3)*(N-4)
assert H == (-8*N**3+(13*m+59)*N*N+(12*m*m-127*m-93)*N
             -5*m**3-18*m*m+245*m-30)
shifted = H.sub(1,2*m+6+t)
assert shifted == ((7*m+3)*(m*m-33*m-64)
                   -(32*m*m+311*m+249)*t-(35*m+85)*t*t-8*t**3)
assert (m*m-33*m-64) == (m-34)*(m+1)-30
assert (m*m-33*m-64).sub(0,35+r) == r*r+37*r+6
assert all(v>0 for v in shifted.sub(2,P(0)).sub(0,35+r).terms.values())
# For N=5m/2+t, m=2+r, every coefficient of H is negative.
H_upper = H.sub(1,F(5,2)*m+t).sub(0,2+r)
assert all(v<0 for v in H_upper.terms.values())
# Stronger family result: F4<1 when m>=4 and N>=5m/2.
G = 4*(N-2)*(N-3)*(N-4)*(N-m-1)-m*(m-1)*A
G_upper = G.sub(1,F(5,2)*m+t).sub(0,4+r)
assert all(v>0 for v in G_upper.terms.values())
Q = N*N-3*(m+1)*N+m*m+5*m+2
assert Q.sub(1,F(5,2)*m+1) == -F(1,4)*m*(m+2)
assert Q.sub(1,2*m+6) == -m*m+5*m+20
# Mutation control: changing the binomial coefficient 30 to 29 must fail.
mutant_A = derived_A - (m-2)*(M-m+1)
try:
    assert A == mutant_A
except AssertionError:
    mutant = 'EXPECTED-FAIL: coefficient 30 changed to 29'
else:
    raise AssertionError('Mutation unexpectedly passed')
print(json.dumps(dict(status='VERIFIED-COMPUTATIONALLY',
                      exact_symbolic_identities=True,
                      H_shift_coefficients_negative=True,
                      G_shift_coefficients_positive=True,
                      negative_control=mutant),indent=2))
