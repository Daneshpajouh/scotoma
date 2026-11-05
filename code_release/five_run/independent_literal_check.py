#!/usr/bin/env python3
"""Compare the new construction with literal four-letter words at N=6.

This is an incidence-checking fixture, not an extension of the theorem's range.
The word list, per-word weights and Hamming distances are built independently.
"""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as F
from itertools import product
from math import lcm
from pathlib import Path
import json
import verify_k3 as checker


def main():
    C, radius = 2, 2
    literal_prefix_masses = {str('0'*j+'1'*(4-j)): F(1, 5) for j in range(5)}
    for text, numerator in [('001X', 2), ('X011', 2), ('011X', 4),
                            ('X001', 4), ('0X00', 4), ('11X1', 4)]:
        for free in '23':
            literal_prefix_masses[text.replace('X', free)] = F(numerator, 10)
    tails = tuple(product('0123', repeat=C))
    weights = {}
    for text, mass in literal_prefix_masses.items():
        tail_weight = radius - int(any(x in text for x in '23'))
        eligible = [t for t in tails if sum(x != '0' for x in t) == tail_weight]
        for tail in eligible:
            weights[tuple(map(int, text))+tuple(map(int, tail))] = mass/len(eligible)
    assert sum(weights.values(), F(0)) == 5
    common_denominator = lcm(*(v.denominator for v in weights.values()))
    integer_weights = [(w, int(v*common_denominator)) for w, v in weights.items()]
    predicted = {}
    literal_maximum = F(0)
    centres = 0
    for centre in product(range(4), repeat=4+C):
        numerator = sum(v for word, v in integer_weights
                        if sum(a != b for a, b in zip(centre, word)) <= radius)
        value = F(numerator, common_denominator)
        prefix = tuple(min(x, 2) for x in centre[:4])
        p = sum(x != 0 for x in centre[4:])
        key = prefix, p
        if key not in predicted:
            predicted[key] = checker.load(C, radius, prefix, p)
        assert value == predicted[key], (centre, value, predicted[key])
        literal_maximum = max(literal_maximum, value)
        centres += 1
    result = {'status': 'VERIFIED-COMPUTATIONALLY', 'N': 4+C, 'm': radius,
              'literal_centres_checked': centres, 'literal_support_points': len(weights),
              'orbit_classes_compared': len(predicted), 'mass': '5',
              'literal_maximum_load': str(literal_maximum), 'all_loads_agree': True}
    root = Path(__file__).resolve().parent
    (root/'independent_literal_result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
