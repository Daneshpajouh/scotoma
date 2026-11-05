#!/usr/bin/env python3
"""C1, with duplicates collapsed. Two of my own runs disagreed and this is why.

pooled_floor.py used distance EXACTLY 1 and found no triangle, matching the appendix's "largest set
off one coordinate is 2". pooled_floor_w.py used distance <= w, which at w=1 includes distance 0, so
the 5 words shared between guides formed cliques of 3 with their neighbours and the max clique read
as 3. The first script printed "duplicates NOT collapsed" as a label and did not act on it.

A pattern that covers a word covers it once, so the object to cover is the set of DISTINCT pooled
words. Collapsing also means pooling starts with a real head start over the per-guide total, because
a shared word is paid for twice there and once here. That head start is reported separately from any
compression, so the two effects are not confused.

Per-guide optimum: sum over guides of ceil(r_g/(w+1)), since each library is an isometric path.
Pooled at w=1: |distinct| - (maximum matching on distance-exactly-1 edges), there being no triangle.
"""
import os, re
from itertools import combinations
from math import ceil

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.join(HERE, '..', 'multiguide_20261017', 'match_changeseq.json')


def library(g):
    ends = [k for k in range(len(g) - 1) if g[k] != g[k + 1]] + [len(g) - 1]
    return [g[:e] + g[e + 1:] for e in ends]


def ham(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


sp = sorted(set(re.findall(r'"([ACGT]{20})"', open(PANEL).read())))
owners = {}
rs = []
for g in sp:
    L = library(g); rs.append(len(L))
    for w in L:
        owners.setdefault(w, set()).add(g)
W = sorted(owners)
V = len(W)
shared = sum(1 for w in W if len(owners[w]) > 1)
print('distinct pooled words %d   words owned by more than one guide %d' % (V, shared))
print('sum of r_g over guides %d, so %d word-slots collapse to %d distinct words'
      % (sum(rs), sum(rs), V))
print()

idx = {w: i for i, w in enumerate(W)}
e1 = [(i, j) for i, j in combinations(range(V), 2) if ham(W[i], W[j]) == 1]
adj = {i: set() for i in range(V)}
for i, j in e1:
    adj[i].add(j); adj[j].add(i)
tri = sum(1 for i, j in e1 if adj[i] & adj[j])
cross = sum(1 for i, j in e1 if not (owners[W[i]] & owners[W[j]]))
print('distance-exactly-1 edges among DISTINCT words: %d   triangles %d   edges with no shared owner %d'
      % (len(e1), tri, cross))
print('  -> largest set agreeing off one coordinate is %s'
      % ('2, matching the appendix' if tri == 0 else 'at least 3'))
print()


def max_matching(n, es):
    g = {v: [] for v in range(n)}
    for a, b in es:
        g[a].append(b); g[b].append(a)
    match = {}

    def aug(v, seen):
        for u in g[v]:
            if u in seen:
                continue
            seen.add(u)
            if u not in match or aug(match[u], seen):
                match[u] = v; match[v] = u
                return True
        return False
    s = 0
    for v in range(n):
        if v not in match and aug(v, set()):
            s += 1
    return s


M = max_matching(V, e1)
pooled = V - M
pg = sum(ceil(r / 2) for r in rs)
print('=== w = 1 ===')
print('  per-guide optimum, sum ceil(r_g/2)        %d' % pg)
print('  pooled, |distinct| - max matching = %d - %d  %d' % (V, M, pooled))
print('  total saving                              %d' % (pg - pooled))
print('    of which from collapsing shared words   %d' % (sum(rs) - V))
print('    of which from cross-guide pairing       %d' % ((pg - pooled) - (sum(rs) - V)))
print()
if (pg - pooled) - (sum(rs) - V) <= 0:
    print('VERDICT: every pattern pooling saves comes from NOT PAYING TWICE for a word two guides')
    print('share. Cross-guide pairing contributes nothing. Pooling is a deduplication, not a')
    print('compression, and C1 dies as a compression question. Reported as a negative.')
else:
    print('VERDICT: cross-guide pairing saves %d patterns beyond deduplication.'
          % ((pg - pooled) - (sum(rs) - V)))
