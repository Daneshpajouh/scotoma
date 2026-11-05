# Row-generation record, $m=8$, $r=11$: the edge is sharp, not a memory limit

The dense orbit LP at $r=11$ is $10{,}240 \times 10{,}240$, which is about 840 MB as float64 and
grows past 2 GB by $r=12$. So the cell was solved by **row generation**, and this file is the record
a referee needs to see that the result is an optimum and not an out-of-memory.

| quantity | value |
|---|---|
| $N$, $m$, $r$ | 19, 8, 11 |
| active prefix length / free length | 10 / 9 |
| support orbits in the LP | 10,240 candidate, 130 in the final dual |
| centre classes in total | **10,240** |
| centre classes in the final active set | **1,927** |
| rounds of generation | **6** |
| wall time | **143 s**, single core |
| violations in the final exact check | **ZERO, over all 10,240 classes** |
| exact dual mass | $19278067/2000000 = 9.6390335$ |
| exact maximum ball mass | $0.9999997$ |
| $\lceil \text{mass} \rceil$ | 10 |
| required for $K=r$ | mass $> 10$ |
| verdict | **not proved; the row ends at $r=10$** |

## Why zero violations makes this an optimum rather than a partial answer

Row generation solves the LP on a **subset** of the centre constraints. Dropping constraints can only
enlarge the feasible region, so the subset optimum is an **upper bound** on the full optimum. The
final dual was then checked against **every** centre class in exact arithmetic and violated none, so
it is feasible for the **full** problem, which makes the subset optimum also a **lower** bound.
Upper and lower coincide, so $9.6390335$ is the full LP optimum.

$\lceil 9.6390335 \rceil = 10 < 11$. **So no fractional dual proves $K=11$ at $m=8$, and the row's
edge of $r=10$ is a property of the problem rather than of our machinery.**

## What it does not say

It does not say $K \ne 11$. The floor may well hold at $r=11$; what is established is that this
instrument cannot reach it. A separate candidate necessary condition for a sub-$r$ cover, if it
holds generally, would forbid failure below $r=13$ at $m=8$, which would mean $r=11$ and $r=12$ hold
and are merely unproved. That condition is derived for a restricted class of centres and is **not**
established here.

Produced by `code_release/walk_floor_rows.py`.
