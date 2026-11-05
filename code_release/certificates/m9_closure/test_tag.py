"""Brute-force checks for tagged shells and dp_tag."""

from __future__ import annotations

import itertools
import subprocess
from pathlib import Path

from shells import pattern
from tagged import Counter, list_shells, shell_size

ROOT = Path(__file__).resolve().parent
BIN = ROOT / "dp_tag"


def hamming(a, b) -> int:
    return sum(x != y for x, y in zip(a, b))


def iter_shell(N: int, m: int, shell):
    g, j, t, u, fl, fr, tagL, tagR = shell
    P = pattern(j, N)
    hasL = j - g - 1 >= 0
    hasR = j + g < N
    farL = list(range(j - g - 1)) if hasL else list(range(max(0, j - g)))
    farR = list(range(j + g + 1, N)) if hasR else list(range(min(N, j + g), N))

    def letters(tag, pat):
        if tag < 0:
            return [None]
        if tag == 0:
            return [pat]
        if tag == 1:
            return [1 - pat]
        return [2, 3]

    for eL in itertools.combinations(farL, t):
        for fL in itertools.combinations(eL, fl):
            fset = set(fL)
            for eR in itertools.combinations(farR, u):
                for fR in itertools.combinations(eR, fr):
                    rset = set(fR)
                    stars = [i for i in list(eL) if i not in fset] + [i for i in list(eR) if i not in rset]
                    for bits in itertools.product((2, 3), repeat=len(stars)):
                        for tl in letters(tagL, 1):
                            for tr in letters(tagR, 0):
                                w = list(P)
                                if tl is not None:
                                    w[j - g - 1] = tl
                                if tr is not None:
                                    w[j + g] = tr
                                for i in fset:
                                    w[i] = 1 - P[i]
                                for i in rset:
                                    w[i] = 1 - P[i]
                                for i, bit in zip(stars, bits):
                                    w[i] = bit
                                yield tuple(w)


def brute_max(N, m, shells, coeffs, alphabet):
    groups = []
    for shell, coeff in zip(shells, coeffs):
        if coeff:
            groups.append((list(iter_shell(N, m, shell)), coeff))
    best = -1
    centre_best = None
    ncent = 0
    for centre in itertools.product(alphabet, repeat=N):
        ncent += 1
        load = 0
        for words, coeff in groups:
            for w in words:
                if hamming(w, centre) <= m:
                    load += coeff
        if load > best:
            best = load
            centre_best = centre
    return best, centre_best, ncent, groups


def write_instance(path: Path, N, m, shells, coeffs):
    rows = []
    for shell, coeff in zip(shells, coeffs):
        if coeff:
            rows.append((*shell, coeff))
    lines = [f"{N} {m} {len(rows)}"]
    for row in rows:
        lines.append(" ".join(str(x) for x in row))
    path.write_text("\n".join(lines) + "\n")


def run_dp(path: Path):
    proc = subprocess.run([str(BIN), str(path)], check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout)
    maximum = None
    witness = None
    classes = None
    for line in proc.stdout.splitlines():
        parts = line.split()
        if parts[0] == "maximum":
            maximum = int(parts[1])
        elif parts[0] == "witness":
            witness = parts[1]
        elif parts[0] == "classes":
            classes = int(parts[1])
    return maximum, witness, classes


def main():
    failures = []
    # Size, anchor distance, and python count against full enumeration.
    N, m = 4, 2
    ctr = Counter(N, m, (1, 2), 2)
    for shell in ctr.shells:
        words = list(iter_shell(N, m, shell))
        if len(words) != shell_size(N, m, *shell):
            failures.append(f"size {shell} enum={len(words)} formula={shell_size(N, m, *shell)}")
            continue
        anchor = tuple(pattern(shell[1], N))
        if any(hamming(w, anchor) != m for w in words):
            failures.append(f"anchor distance {shell}")
        for centre in (anchor, tuple(0 for _ in range(N)), tuple(2 for _ in range(N))):
            got = ctr.count(centre, ctr.shells.index(shell))
            exp = sum(hamming(w, centre) <= m for w in words)
            if got != exp:
                failures.append(f"count {shell} centre={centre} py={got} enum={exp}")
                break
    print(f"enumerated shells={len(ctr.shells)} N={N} m={m}")

    # Every shell at a slightly larger size has anchor count equal to its size.
    big = Counter(6, 3, (1, 2), 2)
    for i, shell in enumerate(big.shells):
        anchor = tuple(pattern(shell[1], 6))
        got = big.count(anchor, i)
        if got != big.sizes[i]:
            failures.append(f"anchor count {shell} got={got} size={big.sizes[i]}")
    print(f"anchor identity shells={len(big.shells)}")

    specs = [
        (4, 2, (1,), 2, 3),
        (5, 2, (1, 2), 2, 4),
        (5, 3, (1, 2), 2, 3),
    ]
    for N, m, guards, flip_cap, nkeep in specs:
        ctr = Counter(N, m, guards, flip_cap)
        shells = ctr.shells[:nkeep]
        coeffs = [3 + i for i in range(len(shells))]
        inst = ROOT / "logs" / f"tag_N{N}_m{m}.in"
        inst.parent.mkdir(exist_ok=True)
        write_instance(inst, N, m, shells, coeffs)
        got, witness, classes = run_dp(inst)
        best, centre, ncent, groups = brute_max(N, m, shells, coeffs, range(3))
        if got != best or classes != 3**N:
            failures.append(f"dp N={N} m={m} dp={got} brute={best} classes={classes}")
        else:
            print(f"PASS dp N={N} m={m} max={got} classes={classes} witness={witness}")
        best4, centre4, n4, _ = brute_max(N, m, shells, coeffs, range(4))
        if best4 != best:
            failures.append(f"symmetry N={N} four={best4} ternary={best} centre={centre4}")
        else:
            print(f"PASS symmetry N={N} m={m} four_classes={n4}")
        if centre4 is not None:
            folded = tuple(2 if x == 3 else x for x in centre4)
            load_f = 0
            load_c = 0
            for words, coeff in groups:
                for w in words:
                    if hamming(w, folded) <= m:
                        load_f += coeff
                    if hamming(w, centre4) <= m:
                        load_c += coeff
            if load_f != load_c:
                failures.append(f"fold mismatch {centre4}")
        # Large coefficients must scale the maximum exactly.
        scale = 10**20
        big_coeffs = [c * scale for c in coeffs]
        write_instance(inst, N, m, shells, big_coeffs)
        got_big, _, _ = run_dp(inst)
        if got_big != best * scale:
            failures.append(f"scale N={N} got={got_big} expected={best * scale}")
        else:
            print(f"PASS scale N={N} m={m}")

    if failures:
        print("FAILURES")
        for item in failures:
            print(item)
        raise SystemExit(1)
    print("ALL PASS")


if __name__ == "__main__":
    main()
