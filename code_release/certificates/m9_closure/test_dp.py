"""Independent checks for the guarded-shell load.

The gold standard enumerates shell words. A centre load is the total weight
of those words at Hamming distance at most m. Letter-swap symmetry is checked
by comparing a centre that uses letter 3 with the centre that rewrites 3 as 2.
"""

from __future__ import annotations

import itertools
import subprocess
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BIN = ROOT / "dp"
M_ALPHA = 4


def pattern(j: int, N: int) -> tuple[int, ...]:
    return tuple(1 if i < j else 0 for i in range(N))


def shells(N: int, m: int, guards: tuple[int, ...], r: int | None = None):
    if r is None:
        r = N + 1
    out = []
    for g in guards:
        for j in range(r):
            left = list(range(max(0, j - g)))
            right = list(range(min(N, j + g), N))
            L, R = len(left), len(right)
            for t in range(m + 1):
                s = m - t
                if t > L or s > R:
                    continue
                for fl in range(min(g - 1, t) + 1):
                    for fr in range(min(g - 1, s) + 1):
                        out.append((g, j, t, fl, fr, left, right))
    return out


def shell_size(N: int, m: int, g: int, j: int, t: int, fl: int, fr: int) -> int:
    L = max(0, j - g)
    R = max(0, N - j - g)
    s = m - t
    if t > L or s > R or fl > t or fr > s or fl > g - 1 or fr > g - 1:
        return 0
    return comb(L, t) * comb(t, fl) * comb(R, s) * comb(s, fr) * (1 << (m - fl - fr))


def enumerate_shell(N: int, m: int, g: int, j: int, t: int, fl: int, fr: int):
    P = pattern(j, N)
    left = list(range(max(0, j - g)))
    right = list(range(min(N, j + g), N))
    s = m - t
    for eL in itertools.combinations(left, t):
        for fL in itertools.combinations(eL, fl):
            fset = set(fL)
            for eR in itertools.combinations(right, s):
                for fR in itertools.combinations(eR, fr):
                    rset = set(fR)
                    stars = [i for i in eL if i not in fset] + [i for i in eR if i not in rset]
                    for bits in itertools.product((2, 3), repeat=len(stars)):
                        w = list(P)
                        for i in fset:
                            w[i] = 1 - P[i]
                        for i in rset:
                            w[i] = 1 - P[i]
                        for i, bit in zip(stars, bits):
                            w[i] = bit
                        yield tuple(w)


def hamming(a, b) -> int:
    return sum(x != y for x, y in zip(a, b))


def brute_loads(N: int, m: int, weighted: list[tuple], alphabet: range):
    """weighted items are (g,j,t,fl,fr,per_word_weight)."""
    words = []
    for g, j, t, fl, fr, wt in weighted:
        if wt == 0:
            continue
        for w in enumerate_shell(N, m, g, j, t, fl, fr):
            words.append((w, wt))
    best = -1
    best_c = None
    count = 0
    for centre in itertools.product(alphabet, repeat=N):
        count += 1
        load = 0
        for w, wt in words:
            if hamming(w, centre) <= m:
                load += wt
        if load > best:
            best = load
            best_c = centre
    return best, best_c, count, words


def write_instance(path: Path, N: int, m: int, rows: list[tuple[int, int, int, int, int, int]]):
    lines = [f"{N} {m} {len(rows)}"]
    for g, j, t, fl, fr, coeff in rows:
        lines.append(f"{g} {j} {t} {fl} {fr} {coeff}")
    path.write_text("\n".join(lines) + "\n")


def run_dp(path: Path) -> tuple[int, str]:
    proc = subprocess.run([str(BIN), str(path)], check=True, capture_output=True, text=True)
    maximum = None
    witness = None
    for line in proc.stdout.splitlines():
        parts = line.split()
        if parts[0] == "maximum":
            maximum = int(parts[1])
        elif parts[0] == "witness":
            witness = parts[1]
    if maximum is None or witness is None:
        raise RuntimeError(proc.stdout + proc.stderr)
    return maximum, witness


def fold3(word):
    return tuple(2 if x == 3 else x for x in word)


def main():
    failures = []
    # Size formula matches enumeration.
    for N, m, guards in ((5, 2, (1, 2)), (6, 3, (1, 2)), (7, 3, (1, 2, 3))):
        for g, j, t, fl, fr, _, _ in shells(N, m, guards):
            got = sum(1 for _ in enumerate_shell(N, m, g, j, t, fl, fr))
            exp = shell_size(N, m, g, j, t, fl, fr)
            if got != exp:
                failures.append(f"size N={N} {(g, j, t, fl, fr)} enum={got} formula={exp}")
    # Every enumerated shell word is at distance m from its anchor and,
    # when the flip budget is at most g-1, at distance at least m+1 from
    # every other full-chain pattern.
    N, m = 6, 3
    for g, j, t, fl, fr, _, _ in shells(N, m, (1, 2)):
        P = [pattern(i, N) for i in range(N + 1)]
        for w in enumerate_shell(N, m, g, j, t, fl, fr):
            if hamming(w, P[j]) != m:
                failures.append(f"anchor distance {(g, j, t, fl, fr)}")
                break
            if any(i != j and hamming(w, P[i]) <= m for i in range(N + 1)):
                failures.append(f"privacy {(g, j, t, fl, fr)} word={w}")
                break
    # Exact DP against full enumeration, including overlapping shells.
    specs = [
        (5, 2, (1,), [(1, 0, 0, 0, 0, 3), (1, 2, 1, 0, 0, 5), (1, 5, 2, 0, 0, 4)]),
        (5, 2, (1, 2), [(2, 0, 0, 0, 1, 2), (1, 2, 1, 0, 0, 7), (2, 0, 0, 0, 1, 1)]),
        (6, 3, (1, 2), [(1, 2, 1, 0, 0, 4), (2, 1, 0, 0, 1, 3), (2, 1, 0, 0, 1, 6)]),
        (6, 3, (2, 3), [(3, 0, 0, 0, 1, 5), (2, 0, 0, 0, 0, 2), (3, 6, 3, 2, 0, 8)]),
    ]
    for N, m, guards, rows in specs:
        present = []
        for g, j, t, fl, fr, coeff in rows:
            if shell_size(N, m, g, j, t, fl, fr) <= 0:
                failures.append(f"illegal test row {(N, m, g, j, t, fl, fr)}")
                continue
            present.append((g, j, t, fl, fr, coeff))
        # Per-word weight equals coeff, so H/S is absorbed and S cancels: the
        # DP coefficient is the per-word weight.
        inst = ROOT / "logs" / f"self_N{N}_m{m}.in"
        inst.parent.mkdir(exist_ok=True)
        write_instance(inst, N, m, present)
        got, witness = run_dp(inst)
        best, centre, ncent, _ = brute_loads(
            N, m, present, range(3)
        )
        if got != best:
            failures.append(f"dp N={N} m={m} dp={got} brute={best} witness={witness} brute_centre={centre}")
        else:
            print(f"PASS dp N={N} m={m} guards={guards} max={got} classes={ncent} witness={witness}")
        # Symmetry: folding letter 3 to 2 preserves the load of every centre.
        best4, centre4, n4, words = brute_loads(N, m, present, range(4))
        if best4 != best:
            failures.append(f"symmetry gap four={best4} ternary={best}")
        else:
            print(f"PASS symmetry N={N} m={m} four_classes={n4} max={best4}")
        # The folded image of the four-letter maximiser matches some ternary load.
        if centre4 is not None:
            folded = fold3(centre4)
            load_f = sum(wt for w, wt in words if hamming(w, folded) <= m)
            load_c = sum(wt for w, wt in words if hamming(w, centre4) <= m)
            if load_f != load_c:
                failures.append(f"fold mismatch {centre4}")
    if failures:
        print("FAILURES")
        for f in failures:
            print(f)
        raise SystemExit(1)
    print("ALL PASS")


if __name__ == "__main__":
    main()
