# THREE ARMS, so that "CFD has skill on these rows" and "CFD's GAP PENALTIES have skill" are
# separated. Each bulge-bearing row carries one gap and several mismatches, so CFD can reach a good
# AUROC through its mismatch terms alone with the 152 gap entries contributing nothing.
#
#   arm B  CFD with the gap column DROPPED, mismatches scored at their true guide positions.
#          This is what a tool that ignores bulges computes.
#   arm C  arm B plus one factor for the gap, taken from the table. For an RNA bulge the unpaired
#          base is a GUIDE base with a definite position, so there is one answer. For a DNA bulge
#          the unpaired base is a TARGET base with no guide position, so it is attributed to a
#          flanking guide position and BOTH conventions are reported; their difference is the
#          measurable cost of the convention.
#
# C minus B is the gap penalties' own contribution. Positions are always computed as the number of
# real GUIDE bases up to the column, so the 24-character DNA-bulge alignments are indexed 1..20 like
# everything else; the earlier run indexed them 1..21 and is withdrawn.
import csv, sys, bisect
from collections import defaultdict
csv.field_size_limit(10 ** 7)
if len(sys.argv) < 2 or sys.argv[1].startswith('-'):
    sys.exit(__doc__ or 'usage: cfd_gap_arms.py CHANGEseq_CR_Lazzarotto_2020_dataset.csv')
P = sys.argv[1]

MM_SCORES = {'rU:dT,12': 0.8, 'rU:dT,13': 0.692307692, 'rU:dC,5': 0.64, 'rG:dA,14': 0.26666666699999997, 'rG:dG,19': 0.448275862, 'rG:dG,18': 0.47619047600000003, 'rG:dG,15': 0.272727273, 'rG:dG,14': 0.428571429, 'rG:dG,17': 0.235294118, 'rG:dG,16': 0.0, 'rC:dC,20': 0.058823529000000006, 'rG:dT,20': 0.9375, 'rG:dG,13': 0.42105263200000004, 'rG:dG,12': 0.529411765, 'rU:dC,6': 0.571428571, 'rU:dG,14': 0.28571428600000004, 'rU:dT,18': 0.666666667, 'rA:dG,13': 0.21052631600000002, 'rA:dG,12': 0.263157895, 'rA:dG,11': 0.4, 'rA:dG,10': 0.333333333, 'rA:dA,19': 0.538461538, 'rA:dA,18': 0.5, 'rA:dG,15': 0.272727273, 'rA:dG,14': 0.214285714, 'rA:dA,15': 0.2, 'rA:dA,14': 0.533333333, 'rA:dA,17': 0.133333333, 'rA:dA,16': 0.0, 'rA:dA,11': 0.307692308, 'rA:dA,10': 0.882352941, 'rA:dA,13': 0.3, 'rA:dA,12': 0.333333333, 'rG:dA,13': 0.3, 'rG:dA,12': 0.384615385, 'rG:dA,11': 0.384615385, 'rG:dA,10': 0.8125, 'rG:dA,17': 0.25, 'rG:dA,16': 0.0, 'rG:dA,15': 0.14285714300000002, 'rG:dA,6': 0.666666667, 'rG:dG,20': 0.428571429, 'rG:dA,19': 0.666666667, 'rG:dA,18': 0.666666667, 'rU:dC,4': 0.625, 'rG:dT,12': 0.933333333, 'rG:dT,13': 0.923076923, 'rU:dG,11': 0.666666667, 'rC:dA,3': 0.6875, 'rC:dA,2': 0.9090909090000001, 'rC:dA,1': 1.0, 'rC:dA,7': 0.8125, 'rC:dA,6': 0.9285714290000001, 'rC:dA,5': 0.636363636, 'rC:dA,4': 0.8, 'rC:dA,9': 0.875, 'rC:dA,8': 0.875, 'rU:dT,6': 0.8666666670000001, 'rA:dG,20': 0.22727272699999998, 'rG:dT,18': 0.692307692, 'rU:dG,10': 0.533333333, 'rG:dT,19': 0.7142857140000001, 'rG:dA,20': 0.7, 'rC:dT,20': 0.5, 'rU:dC,2': 0.84, 'rG:dG,10': 0.4, 'rC:dA,17': 0.46666666700000003, 'rC:dA,16': 0.307692308, 'rC:dA,15': 0.066666667, 'rC:dA,14': 0.7333333329999999, 'rC:dA,13': 0.7, 'rC:dA,12': 0.538461538, 'rC:dA,11': 0.307692308, 'rC:dA,10': 0.9411764709999999, 'rG:dG,11': 0.428571429, 'rU:dC,20': 0.176470588, 'rG:dG,3': 0.384615385, 'rC:dA,19': 0.46153846200000004, 'rC:dA,18': 0.642857143, 'rU:dG,17': 0.705882353, 'rU:dG,16': 0.666666667, 'rU:dG,15': 0.272727273, 'rG:dG,2': 0.692307692, 'rU:dG,13': 0.7894736840000001, 'rU:dG,12': 0.947368421, 'rG:dA,9': 0.533333333, 'rG:dA,8': 0.625, 'rG:dA,7': 0.571428571, 'rG:dG,5': 0.7857142859999999, 'rG:dA,5': 0.3, 'rG:dA,4': 0.363636364, 'rG:dA,3': 0.5, 'rG:dA,2': 0.636363636, 'rG:dA,1': 1.0, 'rG:dG,4': 0.529411765, 'rG:dG,1': 0.7142857140000001, 'rA:dC,9': 0.666666667, 'rG:dG,7': 0.6875, 'rG:dT,5': 0.8666666670000001, 'rU:dT,20': 0.5625, 'rC:dC,15': 0.05, 'rC:dC,14': 0.0, 'rC:dC,17': 0.058823529000000006, 'rC:dC,16': 0.153846154, 'rC:dC,11': 0.25, 'rC:dC,10': 0.38888888899999996, 'rC:dC,13': 0.13636363599999998, 'rC:dC,12': 0.444444444, 'rC:dA,20': 0.3, 'rC:dC,19': 0.125, 'rC:dC,18': 0.133333333, 'rA:dA,1': 1.0, 'rA:dA,3': 0.705882353, 'rA:dA,2': 0.727272727, 'rA:dA,5': 0.363636364, 'rA:dA,4': 0.636363636, 'rA:dA,7': 0.4375, 'rA:dA,6': 0.7142857140000001, 'rA:dA,9': 0.6, 'rA:dA,8': 0.428571429, 'rU:dG,20': 0.090909091, 'rC:dC,9': 0.6190476189999999, 'rC:dC,8': 0.642857143, 'rU:dT,10': 0.857142857, 'rU:dT,11': 0.75, 'rU:dT,16': 0.9090909090000001, 'rU:dT,17': 0.533333333, 'rU:dT,14': 0.6190476189999999, 'rU:dT,15': 0.578947368, 'rC:dC,1': 0.913043478, 'rU:dT,3': 0.7142857140000001, 'rC:dC,3': 0.5, 'rC:dC,2': 0.695652174, 'rC:dC,5': 0.6, 'rC:dC,4': 0.5, 'rC:dC,7': 0.470588235, 'rC:dC,6': 0.5, 'rU:dT,4': 0.47619047600000003, 'rU:dT,8': 0.8, 'rU:dT,9': 0.9285714290000001, 'rA:dC,19': 0.375, 'rA:dC,18': 0.4, 'rA:dC,17': 0.176470588, 'rA:dC,16': 0.192307692, 'rA:dC,15': 0.65, 'rA:dC,14': 0.46666666700000003, 'rA:dC,13': 0.6521739129999999, 'rA:dC,12': 0.7222222220000001, 'rA:dC,11': 0.65, 'rA:dC,10': 0.5555555560000001, 'rU:dC,7': 0.588235294, 'rC:dT,8': 0.65, 'rC:dT,9': 0.857142857, 'rC:dT,6': 0.9285714290000001, 'rC:dT,7': 0.75, 'rC:dT,4': 0.842105263, 'rC:dT,5': 0.571428571, 'rC:dT,2': 0.727272727, 'rC:dT,3': 0.8666666670000001, 'rC:dT,1': 1.0, 'rA:dC,8': 0.7333333329999999, 'rU:dT,1': 1.0, 'rU:dC,3': 0.5, 'rU:dC,1': 0.956521739, 'rU:dT,2': 0.846153846, 'rU:dG,19': 0.275862069, 'rG:dT,14': 0.75, 'rG:dT,15': 0.9411764709999999, 'rG:dT,16': 1.0, 'rG:dT,17': 0.933333333, 'rG:dT,10': 0.933333333, 'rG:dT,11': 1.0, 'rA:dG,9': 0.571428571, 'rA:dG,8': 0.428571429, 'rA:dG,7': 0.4375, 'rA:dG,6': 0.454545455, 'rA:dG,5': 0.5, 'rA:dG,4': 0.352941176, 'rA:dG,3': 0.428571429, 'rA:dG,2': 0.7857142859999999, 'rA:dG,1': 0.857142857, 'rU:dT,5': 0.5, 'rG:dT,2': 0.846153846, 'rA:dC,3': 0.611111111, 'rA:dC,20': 0.764705882, 'rG:dT,1': 0.9, 'rG:dT,6': 1.0, 'rG:dT,7': 1.0, 'rG:dT,4': 0.9, 'rC:dT,19': 0.428571429, 'rG:dG,9': 0.538461538, 'rG:dG,8': 0.615384615, 'rG:dT,8': 1.0, 'rG:dT,9': 0.642857143, 'rU:dG,18': 0.428571429, 'rU:dT,7': 0.875, 'rG:dG,6': 0.681818182, 'rA:dA,20': 0.6, 'rU:dC,9': 0.6190476189999999, 'rA:dG,17': 0.176470588, 'rU:dC,8': 0.7333333329999999, 'rA:dG,16': 0.0, 'rA:dG,19': 0.20689655199999998, 'rG:dT,3': 0.75, 'rU:dG,3': 0.428571429, 'rU:dG,2': 0.857142857, 'rU:dG,1': 0.857142857, 'rA:dG,18': 0.19047619, 'rU:dG,7': 0.6875, 'rU:dG,6': 0.9090909090000001, 'rU:dG,5': 1.0, 'rU:dG,4': 0.647058824, 'rU:dG,9': 0.923076923, 'rU:dG,8': 1.0, 'rU:dC,19': 0.25, 'rU:dC,18': 0.333333333, 'rU:dC,13': 0.260869565, 'rU:dC,12': 0.5, 'rU:dC,11': 0.4, 'rU:dC,10': 0.5, 'rU:dC,17': 0.117647059, 'rU:dC,16': 0.346153846, 'rU:dC,15': 0.05, 'rU:dC,14': 0.0, 'rC:dT,10': 0.8666666670000001, 'rC:dT,11': 0.75, 'rC:dT,12': 0.7142857140000001, 'rC:dT,13': 0.384615385, 'rC:dT,14': 0.35, 'rC:dT,15': 0.222222222, 'rC:dT,16': 1.0, 'rC:dT,17': 0.46666666700000003, 'rC:dT,18': 0.538461538, 'rA:dC,2': 0.8, 'rA:dC,1': 1.0, 'rA:dC,7': 0.705882353, 'rA:dC,6': 0.7142857140000001, 'rA:dC,5': 0.72, 'rA:dC,4': 0.625, 'rU:dT,19': 0.28571428600000004, 'r-:dA,2': 0.692307692308, 'r-:dA,3': 0.4375, 'r-:dA,4': 0.15, 'r-:dA,5': 0.133333333333, 'r-:dA,6': 0.0666666666667, 'r-:dA,7': 0.0, 'r-:dA,8': 0.0, 'r-:dA,9': 0.142857142857, 'r-:dA,10': 0.4, 'r-:dA,11': 0.3125, 'r-:dA,12': 0.0714285714286, 'r-:dA,13': 0.0769230769231, 'r-:dA,14': 0.047619047619, 'r-:dA,15': 0.0526315789474, 'r-:dA,16': 0.0, 'r-:dA,17': 0.0, 'r-:dA,18': 0.153846153846, 'r-:dA,19': 0.142857142857, 'r-:dA,20': 0.6, 'r-:dC,2': 0.714285714286, 'r-:dC,3': 0.461538461538, 'r-:dC,4': 0.176470588235, 'r-:dC,5': 0.0714285714286, 'r-:dC,6': 0.0, 'r-:dC,7': 0.0, 'r-:dC,8': 0.0, 'r-:dC,9': 0.0, 'r-:dC,10': 0.133333333333, 'r-:dC,11': 0.133333333333, 'r-:dC,12': 0.0, 'r-:dC,13': 0.0, 'r-:dC,14': 0.0, 'r-:dC,15': 0.0, 'r-:dC,16': 0.0, 'r-:dC,17': 0.0588235294118, 'r-:dC,18': 0.142857142857, 'r-:dC,19': 0.0357142857143, 'r-:dC,20': 0.0454545454545, 'r-:dT,2': 0.727272727273, 'r-:dT,3': 0.6, 'r-:dT,4': 0.1, 'r-:dT,5': 0.0, 'r-:dT,6': 0.0833333333333, 'r-:dT,7': 0.0, 'r-:dT,8': 0.0625, 'r-:dT,9': 0.1875, 'r-:dT,10': 0.25, 'r-:dT,11': 0.0769230769231, 'r-:dT,12': 0.0769230769231, 'r-:dT,13': 0.0, 'r-:dT,14': 0.0, 'r-:dT,15': 0.0, 'r-:dT,16': 0.0, 'r-:dT,17': 0.0, 'r-:dT,18': 0.0, 'r-:dT,19': 0.307692307692, 'r-:dT,20': 0.6, 'r-:dG,2': 0.96, 'r-:dG,3': 0.666666666667, 'r-:dG,4': 0.375, 'r-:dG,5': 0.24, 'r-:dG,6': 0.142857142857, 'r-:dG,7': 0.0588235294118, 'r-:dG,8': 0.133333333333, 'r-:dG,9': 0.238095238095, 'r-:dG,10': 0.333333333333, 'r-:dG,11': 0.1, 'r-:dG,12': 0.0555555555556, 'r-:dG,13': 0.0434782608696, 'r-:dG,14': 0.0, 'r-:dG,15': 0.0, 'r-:dG,16': 0.0, 'r-:dG,17': 0.0, 'r-:dG,18': 0.2, 'r-:dG,19': 0.4375, 'r-:dG,20': 0.529411764706, 'rA:d-,2': 0.969230769231, 'rA:d-,3': 0.888888888889, 'rA:d-,4': 0.58064516129, 'rA:d-,5': 0.328125, 'rA:d-,6': 0.276923076923, 'rA:d-,7': 0.430769230769, 'rA:d-,8': 0.234375, 'rA:d-,9': 0.276923076923, 'rA:d-,10': 0.569230769231, 'rA:d-,11': 0.276923076923, 'rA:d-,12': 0.0615384615385, 'rA:d-,13': 0.03125, 'rA:d-,14': 0.0461538461538, 'rA:d-,15': 0.03125, 'rA:d-,16': 0.0153846153846, 'rA:d-,17': 0.0, 'rA:d-,18': 0.0, 'rA:d-,19': 0.0, 'rA:d-,20': 0.0615384615385, 'rC:d-,2': 0.96875, 'rC:d-,3': 0.933333333333, 'rC:d-,4': 0.68253968254, 'rC:d-,5': 0.365079365079, 'rC:d-,6': 0.4375, 'rC:d-,7': 0.569230769231, 'rC:d-,8': 0.123076923077, 'rC:d-,9': 0.415384615385, 'rC:d-,10': 0.615384615385, 'rC:d-,11': 0.169230769231, 'rC:d-,12': 0.046875, 'rC:d-,13': 0.0, 'rC:d-,14': 0.0461538461538, 'rC:d-,15': 0.0, 'rC:d-,16': 0.0, 'rC:d-,17': 0.0, 'rC:d-,18': 0.0, 'rC:d-,19': 0.0, 'rC:d-,20': 0.03125, 'rU:d-,2': 0.953125, 'rU:d-,3': 0.921875, 'rU:d-,4': 0.65625, 'rU:d-,5': 0.215384615385, 'rU:d-,6': 0.446153846154, 'rU:d-,7': 0.553846153846, 'rU:d-,8': 0.296875, 'rU:d-,9': 0.609375, 'rU:d-,10': 0.650793650794, 'rU:d-,11': 0.323076923077, 'rU:d-,12': 0.0307692307692, 'rU:d-,13': 0.0, 'rU:d-,14': 0.0615384615385, 'rU:d-,15': 0.0, 'rU:d-,16': 0.0, 'rU:d-,17': 0.0, 'rU:d-,18': 0.0, 'rU:d-,19': 0.0, 'rU:d-,20': 0.0, 'rG:d-,2': 0.921875, 'rG:d-,3': 0.822580645161, 'rG:d-,4': 0.523076923077, 'rG:d-,5': 0.296875, 'rG:d-,6': 0.390625, 'rG:d-,7': 0.492063492063, 'rG:d-,8': 0.301587301587, 'rG:d-,9': 0.328125, 'rG:d-,10': 0.587301587302, 'rG:d-,11': 0.25, 'rG:d-,12': 0.125, 'rG:d-,13': 0.046875, 'rG:d-,14': 0.0634920634921, 'rG:d-,15': 0.0, 'rG:d-,16': 0.0, 'rG:d-,17': 0.0, 'rG:d-,18': 0.0, 'rG:d-,19': 0.0634920634921, 'rG:d-,20': 0.09375}
PAM_SCORES = {'AA': 0.0, 'AC': 0.0, 'GT': 0.016129031999999998, 'AG': 0.25925925899999996, 'CC': 0.0, 'CA': 0.0, 'CG': 0.107142857, 'TT': 0.0, 'GG': 1.0, 'GC': 0.022222222000000003, 'AT': 0.0, 'GA': 0.06944444400000001, 'TG': 0.038961038999999996, 'TA': 0.0, 'TC': 0.0, 'CT': 0.0}

MM, PM = MM_SCORES, PAM_SCORES
RC = {'A':'T','C':'G','G':'C','T':'A','U':'A','-':'-'}

import random as _rnd
_r = _rnd.Random(20261022)
_shuf = list(range(1, 21)); _r.shuffle(_shuf)
PERM = {p: _shuf[p-1] for p in range(1, 21)}
print('permutation seed 20261022 ->', [PERM[p] for p in range(1,21)], flush=True)


def arms(g_al, t_al, pam):
    """(score_B, score_C_after, score_C_before, ok). g_al/t_al are spacer-only, equal length."""
    gpos = 0          # real guide bases seen so far
    base_mm = 1.0
    gap_pos_after = gap_pos_before = None
    gap_dna = gap_guide = None
    for i, (gc, tc) in enumerate(zip(g_al, t_al)):
        if gc == '-':                      # DNA bulge: target base unpaired, no guide position
            if tc not in RC: return 0.0, 0.0, 0.0, 0.0, 0.0, 'iupac'
            gap_dna = tc
            gap_pos_before = gpos          # position of the guide base before the gap
            gap_pos_after = gpos + 1       # position of the guide base after it
            continue
        gpos += 1
        if tc == '-':                      # RNA bulge: THIS guide base is unpaired
            gap_guide = gc
            gap_pos_after = gap_pos_before = gpos
            continue
        if gc == tc:
            continue
        if tc not in RC or gc not in RC:
            # CRISPRme's own behaviour: revcom returns None on an untranslatable base, the key
            # concatenation raises, and the handler sets score 0 and breaks. An N in the target
            # therefore zeroes the row rather than skipping it. Matched here, and counted.
            return 0.0, 0.0, 0.0, 0.0, 0.0, 'iupac'
        k = 'r' + gc.replace('T','U') + ':d' + RC[tc] + ',' + str(gpos)
        v = MM.get(k)
        if v is None:
            # INVARIANT, not a diagnostic. Every lookup must hit the table or this is not CFD.
            # Previously this returned a flag that was counted and printed and read past.
            raise AssertionError('CFD table miss at %r: the scoring is not CFD' % k)
        base_mm *= v
    pam_f = PM.get(pam, 0.0)
    B = base_mm * pam_f
    def gap_factor(p):
        if p is None or not (1 <= p <= 20):
            return None
        if gap_guide is not None:                    # rX:d-
            return MM.get('r' + gap_guide.replace('T','U') + ':d-,' + str(p))
        if gap_dna is not None:                      # r-:dX
            return MM.get('r-:d' + RC[gap_dna] + ',' + str(p))
        return 1.0
    fa, fb = gap_factor(gap_pos_after), gap_factor(gap_pos_before)
    # R and PI reindex the SAME penalty table. The gap's identity and orientation are unchanged;
    # only the position the penalty is looked up at moves. Reversal is 21 - p; the permutation is
    # a fixed shuffle of 1..20 so the arm is reproducible from the seed alone.
    pa = gap_pos_after
    fr = gap_factor(21 - pa) if pa is not None and 1 <= 21 - pa <= 20 else None
    fp = gap_factor(PERM[pa]) if pa is not None and 1 <= pa <= 20 else None
    return (B,
            (B * fa if fa is not None else None),
            (B * fb if fb is not None else None),
            (B * fr if fr is not None else None),
            (B * fp if fp is not None else None),
            True)

def auroc(pairs):
    pos = [s for s, l in pairs if l == 1]
    neg = sorted(s for s, l in pairs if l == 0)
    if not pos or not neg: return None
    tot = 0.0
    for p in pos:
        lo = bisect.bisect_left(neg, p); hi = bisect.bisect_right(neg, p)
        tot += lo + 0.5 * (hi - lo)
    return tot / (len(pos) * len(neg))

strata = defaultdict(list)
bad = scored = iupac = 0
with open(P, newline='') as fh:
    for row in csv.DictReader(fh):
        g = (row.get('Align.sgRNA') or '').strip().upper()
        t = (row.get('Align.off-target') or '').strip().upper()
        if len(g) != len(t) or len(g) < 23: continue
        try:
            nb = int(float((row.get('Align.#Bulges') or '0').strip()))
            nm = int(float((row.get('Align.#Mismatches') or '0').strip()))
            reads = float((row.get('reads') or '0') or 0)
        except Exception: continue
        if nb != 1: continue                       # single-bulge rows only, which is the class
        if (g.count('-') + t.count('-')) != 1: continue
        B, Ca, Cb, Cr, Cp, ok = arms(g[:-3], t[:-3], t[-2:])
        if ok == 'iupac': iupac += 1
        elif not ok or B is None: bad += 1; continue
        scored += 1
        kind = 'RNA' if t.count('-') else 'DNA'
        strata[(kind, nm + nb)].append((B, Ca, Cb, Cr, Cp, 1 if reads > 0 else 0))

print('single-bulge rows scored: %d   unscorable: %d   zeroed by an IUPAC base: %d'
      % (scored, bad, iupac), flush=True)

# Paired DeLong for the C-versus-B comparison, as pre-registered: the two scores are on the SAME
# samples, so the variance of the difference needs the covariance term and an unpaired interval
# would be wrong. Holm across the strata that clear the floor.
import math

def delong_pair(pos_a, neg_a, pos_b, neg_b):
    """AUCs and the variance of their difference, DeLong 1988, placement-value form."""
    def placements(pos, neg):
        ns = sorted(neg); m, n = len(pos), len(neg)
        v10 = []
        for p in pos:
            lo = bisect.bisect_left(ns, p); hi = bisect.bisect_right(ns, p)
            v10.append((lo + 0.5 * (hi - lo)) / n)
        ps = sorted(pos)
        v01 = []
        for q in neg:
            lo = bisect.bisect_left(ps, q); hi = bisect.bisect_right(ps, q)
            v01.append((m - hi + 0.5 * (hi - lo)) / m)
        return v10, v01
    a10, a01 = placements(pos_a, neg_a)
    b10, b01 = placements(pos_b, neg_b)
    m, n = len(pos_a), len(neg_a)
    Aa, Ab = sum(a10) / m, sum(b10) / m
    def cov(x, y, mean_x, mean_y, k):
        return sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y)) / (k - 1)
    s10_aa = cov(a10, a10, Aa, Aa, m); s10_bb = cov(b10, b10, Ab, Ab, m)
    s10_ab = cov(a10, b10, Aa, Ab, m)
    s01_aa = cov(a01, a01, Aa, Aa, n); s01_bb = cov(b01, b01, Ab, Ab, n)
    s01_ab = cov(a01, b01, Aa, Ab, n)
    var = ((s10_aa + s10_bb - 2 * s10_ab) / m) + ((s01_aa + s01_bb - 2 * s01_ab) / n)
    return Aa, Ab, var

def norm_sf(z):
    return 0.5 * math.erfc(abs(z) / math.sqrt(2)) * 2

print()
print('%-4s %5s %7s %9s  %7s %7s %9s %19s %9s'
      % ('kind','edits','pos','neg','B drop','arm','arm - B','95% CI','raw p'))
tests = []
for key in sorted(strata, key=lambda k: (k[0], k[1])):
    rows = strata[key]
    pos = [r for r in rows if r[5] == 1]; neg = [r for r in rows if r[5] == 0]
    if len(pos) < 20 or len(neg) < 20:
        print('%-4s %5d %7d %9d  (below the pre-registered 20/20 floor)' % (key[0], key[1], len(pos), len(neg)))
        continue
    for lbl, idx in (('C after', 1), ('C before', 2), ('R reversed', 3), ('PI permuted', 4)):
        if key[0] == 'RNA' and lbl == 'C before':
            continue                       # convention-free for RNA; identical to "after"
        # A convention can leave the gap OUTSIDE guide positions 1..20, in which case the site has
        # no CFD value under that convention at all. That is informative rather than a nuisance, so
        # the pair is restricted to rows where BOTH arms are defined and the loss is reported.
        posd = [r for r in pos if r[idx] is not None]
        negd = [r for r in neg if r[idx] is not None]
        lost = (len(pos) - len(posd)) + (len(neg) - len(negd))
        if len(posd) < 20 or len(negd) < 20:
            print('%-4s %5d %7d %9d  the %s convention leaves %d row(s) with NO CFD value, and the '
                  'defined remainder is below the floor' % (key[0], key[1], len(posd), len(negd), lbl, lost))
            continue
        pa = [r[idx] for r in posd]; na = [r[idx] for r in negd]
        pb = [r[0] for r in posd];   nb = [r[0] for r in negd]
        # INVARIANT: both arms must actually vary inside the stratum, or their AUROC is 0.5 by
        # construction and the comparison is a tautology. The first design failed exactly here.
        if len(set(pa + na)) < 2 or len(set(pb + nb)) < 2:
            raise AssertionError('an arm is constant within stratum %r: AUROC would be 0.5 by '
                                 'construction and the comparison is vacuous' % (key,))
        Ac, Ab_, var = delong_pair(pa, na, pb, nb)
        d = Ac - Ab_
        se = math.sqrt(var) if var > 0 else float('nan')
        z = d / se if se == se and se > 0 else float('nan')
        p = norm_sf(z) if z == z else float('nan')
        tests.append((key, lbl, d, p))
        tag = lbl
        print('%-4s %5d %7d %9d  %7.4f %7.4f %+9.4f  [%+.4f, %+.4f] %9.2e  %s%s'
              % (key[0], key[1], len(posd), len(negd), Ab_, Ac, d, d - 1.96 * se, d + 1.96 * se, p, tag,
                 ('  (%d row(s) undefined under this convention)' % lost) if lost else ''))
print()
k = len(tests)
order = sorted(range(k), key=lambda i: tests[i][3])
print('Holm across the %d tests, in p order:' % k)
holm_ok = True
for rank, i in enumerate(order):
    key, lbl, d, p = tests[i]
    thresh = 0.05 / (k - rank)
    mark = 'significant' if p < thresh else 'NOT significant'
    print('   %-4s %d %-12s d=%+0.4f  p=%.3e  vs %.3e  %s' % (key[0], key[1], lbl, d, p, thresh, mark))
