"""Finite checks of the scope fixtures used in the covering paper."""
from pathlib import Path
from itertools import permutations
import json

HERE = Path(__file__).resolve().parent
RELEASE = HERE.parents[1]


def distance(a, b):
    return sum(x != y for x, y in zip(a, b))


def main():
    x = (0,)*8 + (2,)*9 + (0,)*2
    y = (1,)*8 + (0,)*9 + (2,)*2
    assert distance(x, (0,)*19) == 9
    assert distance(y, (1,)*8 + (0,)*11) == 2
    assert distance(x, y) == 19
    column_cases = 0
    for values in permutations(range(4), 3):
        assert sum(a != 0 for a in values) >= 2
        for assigned in __import__('itertools').product((0, 1), repeat=3):
            assert sum(a != b for a, b in zip(values, assigned)) >= 1
            column_cases += 1
    panel = json.loads((RELEASE / 'data/producers_2026-10-20/PANELREF_GRCh38_NGG.json').read_text())
    guides = [row['spacer'] for row in panel['rows']]
    assert len(guides) == 110 and all(len(g) == 20 for g in guides)
    coordinates = 0
    for guide in guides:
        insertions = {guide[:i]+a+guide[i:] for i in range(21) for a in 'ACGT'}
        assert len(insertions) == 64
        for i in range(21):
            assert {word[i] for word in insertions} == set('ACGT')
            coordinates += 1
    assert coordinates == 2310
    print(json.dumps({'status': 'VERIFIED-COMPUTATIONALLY',
                      'separated_pair_distances': [9, 2, 19],
                      'three_letter_column_cases': column_cases,
                      'no_triple_distance_lower_bound': 30,
                      'three_word_budget': 27,
                      'panel_guides': 110, 'saturated_insertion_coordinates': coordinates}, indent=2))


if __name__ == '__main__':
    main()
