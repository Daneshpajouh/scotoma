"""Literal distance checks for the elementary construction and its printed tables."""
import argparse
import json
from pathlib import Path

OFFSETS = ((4, 2), (6, 3), (9, 2), (12, 3), (16, 2))
EXPECTED_TABLES = {
    4: list(zip((9,9,9,11,10,10,13,12,13,14,14,14,14,14,16,16,17,18,15),
                (5,6,7,7,8,8,9,10,10,11,12,12,13,14,15,16,17,18,19))),
    5: list(zip((11,11,11,11,11,11,11,13,14,14,15,12,15,15,17,16,17,18,14),
                (6,7,8,8,9,9,10,11,11,12,13,13,14,15,16,16,17,18,19))),
}


def hamming(a, b):
    assert len(a) == len(b)
    return sum(x != y for x, y in zip(a, b))


def words(width, count, offsets, cycle):
    patterns, witnesses, supports = [], [], []
    for k in range(count):
        p = [1] * k + [0] * (width - k)
        w = p.copy()
        s = {}
        for off, letter in offsets:
            position = (k + off - 1) % cycle
            w[position] = letter
            s[position] = letter
        patterns.append(p)
        witnesses.append(w)
        supports.append(s)
    return patterns, witnesses, supports


def check(width, budget, offsets, cycle):
    p, w, s = words(width, width+1, offsets, cycle)
    radii = {hamming(a, b) for a, b in zip(p, w)}
    pair_min = private_min = width+1
    pairs = directed = 0
    for k in range(width+1):
        for j in range(width+1):
            if k == j:
                continue
            boundary = set(range(min(k, j), max(k, j)))
            actual = hamming(w[k], p[j])
            assert actual == len(s[k]) + len(boundary - s[k].keys())
            private_min = min(private_min, actual)
            directed += 1
            if j > k:
                actual = hamming(w[k], w[j])
                union = s[k].keys() | s[j].keys()
                formula = (len(s[k].keys() ^ s[j].keys())
                           + sum(s[k][q] != s[j][q] for q in s[k].keys() & s[j].keys())
                           + len(boundary - union))
                assert actual == formula
                pair_min = min(pair_min, actual)
                pairs += 1
    valid = radii == {budget} and pair_min > 2*budget and private_min > budget
    return dict(width=width, budget=budget, cycle=cycle, radii=sorted(radii),
                pair_min=pair_min, private_min=private_min, pairs=pairs,
                directed_private=directed, valid=valid)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--appendix', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    results = []
    for m in range(6):
        for n in (20, 21, 25, 30, 40, 60):
            row = check(n-1, m, OFFSETS[:m], 19)
            assert row['valid'], row
            results.append(row)
        print('checked all six lengths, budget', m, flush=True)
    for width in (18, 17):
        row = check(width, 4, ((2, 2), (5, 3), (9, 2), (13, 3)), width)
        assert row['valid'], row
        results.append(row)
        print('checked shortened burst width', width, flush=True)
    table = []
    for m in (4, 5):
        p, w, _ = words(19, 20, OFFSETS[:m], 19)
        values = [(hamming(w[0], w[t]), hamming(w[0], p[t])) for t in range(1,20)]
        assert values == EXPECTED_TABLES[m], (m, values)
        table.append(values)
    print('all 76 recorded table distances match literal words', flush=True)
    if args.appendix:
        section = args.appendix.read_text().split('The nineteen cyclic cases, at both budgets.}', 1)[1]
        rows = []
        for line in section.splitlines():
            cells = [s.strip().replace('\\\\','') for s in line.split('&')]
            if len(cells) == 5 and cells[0].isdigit():
                rows.append([int(s) for s in cells])
        assert len(rows) == 19
        for t, row in enumerate(rows, 1):
            assert row == [t, *table[0][t-1], *table[1][t-1]], (t, row)
        print('all 76 printed table distances match literal words', flush=True)
    bad_letter = list(OFFSETS)
    bad_letter[-1] = (16, 3)
    controls = [check(19, 5, bad_letter, 19), check(19, 5, OFFSETS[1:], 19)]
    assert all(not row['valid'] for row in controls)
    print('rejected both corrupted constructions', flush=True)
    report = dict(status='PASS', checks=results, tables=table,
                  rejected_controls=controls, proof_scope='finite implementation checks; analytic reduction supplies all lengths')
    if args.out:
        args.out.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
