#!/usr/bin/env python3
"""Count the private words of the brute-force check at appendix_moved.tex:1767-1769.

The paper's convention (appendix_moved.tex:1699, "n-1=19" at n = 20): n patterns P_j = 1^j 0^(L-j),
j = 0..n-1, over words of length L = n-1 on four letters {0,1,2,3}. A word x is private to P_j at
distance m when d(x, P_j) == m and d(x, P_k) > m for every other k (Hamming distance).
Cases (n, m) are the five the sentence lists. Every word of length L is enumerated; nothing is sampled.

Usage: python3 -I private_words.py
Prints the count per case and the total. Exit 0. No paper value is typed here.
"""
import itertools

CASES = [(7, 2), (8, 2), (7, 3), (8, 3), (9, 2)]


def dist(a, b):
    return sum(x != y for x, y in zip(a, b))


def private_count(n, m):
    L = n - 1
    pats = [tuple([1] * j + [0] * (L - j)) for j in range(n)]
    count = 0
    for x in itertools.product(range(4), repeat=L):
        near = [p for p in pats if dist(x, p) <= m]
        if len(near) == 1 and dist(x, near[0]) == m:
            count += 1
    return count


total = 0
for n, m in CASES:
    c = private_count(n, m)
    total += c
    print(f"(n, m) = ({n}, {m}): {c} private words at distance m")
print(f"private words across the five cases: {total}")
