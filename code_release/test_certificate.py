#!/usr/bin/env python3
"""Exact numerical check of Theorem 2(i) on random discrete classes. Run: python3 test_certificate.py

The certificate: take the k*-th largest of m1 positive scores as the threshold t_hat and certify
if and only if no negative score reaches t_hat, with k* the least k such that
P(Bin(m1, r) >= k) <= delta2, and m0 the least integer with (1 - gamma)^m0 <= delta1.

Everything is computed exactly with rational arithmetic, with no simulation:
  * P(t_hat >= t) = P(Bin(m1, G1(t)) >= k*), where G1(t) = P_1([t, inf)); this also holds with ties.
  * P(certify) = sum_t P(t_hat = t) * (1 - G0(t))^m0, G0(t) = P_0([t, inf)), by independence.
Checked: (a) the tolerance statement P(G1(t_hat) < r) <= delta2 for random positive laws;
(b) the type-I error is at most delta1 + delta2 for every non-attaining law of a random class,
with gamma computed exactly from that class; (c) the power at attaining laws with rho* > r is at
least P(Bin(m1, rho*) >= k*), and it is close to one for large m1.
"""
import random
import sys
from fractions import Fraction as F
from math import comb, log

from m0_search import m0_iid

SCORES = list(range(8))


def binom_tail_ge(m, p, k):
    """P(Bin(m, p) >= k), exact, for Fraction p."""
    if k <= 0:
        return F(1)
    return sum(F(comb(m, j)) * p ** j * (1 - p) ** (m - j) for j in range(k, m + 1))


def kstar(m1, r, delta2):
    k = 0
    while binom_tail_ge(m1, r, k) > delta2:
        k += 1
        assert k <= m1, "m1 too small for this r and delta2"
    return k


def tails(pmf):
    """G(t) = P([t, inf)) over the integer scores 0..7, plus G(8) = 0."""
    out, acc = {8: F(0)}, F(0)
    for t in reversed(SCORES):
        acc += pmf[t]
        out[t] = acc
    return out


def rho_star(p0, p1):
    """Recall of the best zero-false-alarm upper set, discrete scores. Upper sets are [t, inf)
    and (t, inf) = [t+1, inf); the largest with zero negative mass has recall G1 of its start."""
    g0, g1 = tails(p0), tails(p1)
    best = F(0)
    for t in SCORES + [8]:
        if g0[t] == 0:
            best = max(best, g1[t])
    return best


def gap_of_law(p0, p1, r):
    """inf{ alpha(U) : U upper, beta(U) >= r } for one law, upper sets [t, inf), t = 0..8."""
    g0, g1 = tails(p0), tails(p1)
    vals = [g0[t] for t in SCORES + [8] if g1[t] >= r]
    return min(vals) if vals else None


def random_pmf(rng, support):
    w = [F(rng.randint(0, 6)) if s in support else F(0) for s in SCORES]
    if sum(w) == 0:
        w[support[0]] = F(1)
    z = sum(w)
    return [x / z for x in w]


def prob_certify(p0, p1, m0, m1, ks):
    g0, g1 = tails(p0), tails(p1)
    total = F(0)
    for t in SCORES:
        p_ge_t = binom_tail_ge(m1, g1[t], ks)
        p_ge_t1 = binom_tail_ge(m1, g1[t + 1], ks)
        total += (p_ge_t - p_ge_t1) * (1 - g0[t]) ** m0
    return total          # k* <= m1 so t_hat is always defined


def test_tolerance_statement():
    rng = random.Random(1)
    r, d2, m1 = F(9, 10), F(1, 40), 60
    ks = kstar(m1, r, d2)
    worst, n = F(0), 0
    for _ in range(200):
        p1 = random_pmf(rng, SCORES)
        g1 = tails(p1)
        # P(G1(t_hat) < r) = sum over t with G1(t) < r of P(t_hat = t)
        bad = sum(binom_tail_ge(m1, g1[t], ks) - binom_tail_ge(m1, g1[t + 1], ks)
                  for t in SCORES if g1[t] < r)
        assert bad <= d2, (bad, d2)
        worst, n = max(worst, bad), n + 1
    return f"{n} random positive laws, worst P(G1(t_hat)<r) = {float(worst):.4f} <= {float(d2)}"


def test_level_and_power():
    rng = random.Random(2)
    r, d1, d2, m1 = F(9, 10), F(1, 40), F(1, 40), 60
    ks = kstar(m1, r, d2)
    for trial in range(3):
        laws = []
        while len(laws) < 150:
            p0 = random_pmf(rng, SCORES)
            p1 = random_pmf(rng, SCORES)
            laws.append((p0, p1))
        non_att = [(p0, p1) for p0, p1 in laws if rho_star(p0, p1) < r]
        gaps = [gap_of_law(p0, p1, r) for p0, p1 in non_att]
        gamma = min(g for g in gaps if g is not None)
        assert gamma > 0, "a random class should have a positive gap"
        m0 = m0_iid(gamma, d1)
        worst = max(prob_certify(p0, p1, m0, m1, ks) for p0, p1 in non_att)
        assert worst <= d1 + d2, (trial, float(worst))
    # attaining laws with rho* > r: negatives below, positives above
    for _ in range(40):
        p0 = random_pmf(rng, [0, 1, 2, 3])
        p1 = random_pmf(rng, [4, 5, 6, 7])
        rho = rho_star(p0, p1)
        assert rho > r
        for big in (60, 400):
            kb = kstar(big, r, d2)
            power = prob_certify(p0, p1, m0, big, kb)
            assert power >= binom_tail_ge(big, rho, kb), (float(power),)
    return (f"3 random classes ({len(non_att)} non-attaining laws in the last), gamma={float(gamma):.4f}, "
            f"m0={m0}, worst type-I = {float(worst):.4f} <= {float(d1 + d2)}; "
            f"power bound verified at m1 = 60 and 400")


def test_adversarial_boundary_family():
    """Non-attaining laws with recall just below r and false-alarm mass exactly gamma on the rule
    that meets the recall target: the tightest case for the level bound."""
    r, d1, d2, m1 = F(9, 10), F(1, 40), F(1, 40), 60
    ks = kstar(m1, r, d2)
    out = []
    for gamma in (F(1, 20), F(1, 100)):
        m0 = m0_iid(gamma, d1)
        worst = F(0)
        for q in (F(1, 2), F(4, 5), F(85, 100), F(89, 100), F(899, 1000), F(8999, 10000)):
            # positives: mass q at score 7, rest at score 3; negatives: gamma at score 3, rest at 0
            p1 = [F(0)] * 8
            p1[7], p1[3] = q, 1 - q
            p0 = [F(0)] * 8
            p0[3], p0[0] = gamma, 1 - gamma
            assert rho_star(p0, p1) == q < r and gap_of_law(p0, p1, r) == gamma
            worst = max(worst, prob_certify(p0, p1, m0, m1, ks))
        assert worst <= d1 + d2, (float(gamma), float(worst))
        out.append(f"gamma={float(gamma)}: m0={m0}, worst type-I {float(worst):.4f}")
    return "; ".join(out) + f" (bound {float(d1 + d2)})"


TESTS = [test_tolerance_statement, test_level_and_power, test_adversarial_boundary_family]

if __name__ == "__main__":
    failed = 0
    for t in TESTS:
        try:
            print(f"PASS  {t.__name__}: {t()}")
        except Exception as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(TESTS) - failed} of {len(TESTS)} tests passed")
    sys.exit(1 if failed else 0)
