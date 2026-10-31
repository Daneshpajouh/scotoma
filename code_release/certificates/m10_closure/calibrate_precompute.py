"""Compare grouped and reference precomputation with literal four-letter maxima."""
from itertools import combinations, product
from pathlib import Path
import json
import os
import random
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from support import layout, list_shells, pattern, shell_size


def expand(N, depth, shell):
    j, guard, left, right, t, u, fl, fr = shell
    lay = layout(j, N, depth)
    base = pattern(j, N)
    fixed = list(zip(lay['guard'], guard)) + list(zip(lay['resL'], left)) + list(zip(lay['resR'], right))
    def arm(positions, edits, flips):
        for changed in combinations(positions, edits):
            for flipped in combinations(changed, flips):
                flipped = set(flipped)
                starred = [x for x in changed if x not in flipped]
                for letters in product((2, 3), repeat=len(starred)):
                    yield [(x, 1-base[x]) for x in flipped] + list(zip(starred, letters))
    words = []
    for chosen in product(*[(base[pos],) if symbol == 0 else (1-base[pos],) if symbol == 1 else (2, 3)
                            for pos, symbol in fixed]):
        for a in arm(lay['farL'], t, fl):
            for b in arm(lay['farR'], u, fr):
                word = list(base)
                for pos, value in list(zip([x for x, _ in fixed], chosen)) + a + b:
                    word[pos] = value
                words.append(tuple(word))
    assert len(words) == len(set(words)) == shell_size(N, t+u+sum(x != 0 for _, x in fixed), depth, shell)
    return words


def run(binary, work):
    cases = [(4,0,0,0,False), (4,1,0,1,False), (5,2,1,2,False),
             (5,2,2,1,False), (5,2,2,2,True), (8,1,2,1,False),
             (8,2,2,2,False), (8,2,1,2,True), (8,2,0,2,False)]
    reports = []
    for index, (N,m,depth,cap,free) in enumerate(cases):
        shells = list_shells(N,m,depth,cap,free_guard=free)
        rng = random.Random(9100+index)
        chosen = rng.sample(shells, min(8,len(shells)))
        coefficients = [rng.randrange(1,18) for _ in chosen]
        points = {}
        for shell, coeff in zip(chosen, coefficients):
            for word in expand(N,depth,shell):
                value = sum(x << (2*k) for k,x in enumerate(word))
                points[value] = points.get(value,0)+coeff
        mask = sum(1 << (2*k) for k in range(N))
        def score(centre):
            total = 0
            for word, weight in points.items():
                delta = centre ^ word
                if ((delta | (delta >> 1)) & mask).bit_count() <= m:
                    total += weight
            return total
        literal = max(score(centre) for centre in range(4**N))
        lines = [f'{N} {m} {depth} {cap} {int(free)} {len(chosen)}']
        for shell, coeff in zip(chosen, coefficients):
            j,g,left,right,t,u,fl,fr=shell
            lines.append(' '.join(map(str,[j,len(g),len(left),len(right),*g,*left,*right,t,u,fl,fr,coeff])))
        path = Path(work) / f'calibration_{index}.in'
        path.write_text('\n'.join(lines)+'\n')
        outputs = []
        for reference in (False,True):
            env = os.environ.copy()
            env['DP_THREADS'] = '2'
            env.pop('DP_REFERENCE_PRECOMPUTE',None)
            if reference:
                env['DP_REFERENCE_PRECOMPUTE'] = '1'
            result = subprocess.run([str(binary),'max',str(path)], env=env, text=True, capture_output=True, check=True)
            fields = dict(line.split(maxsplit=1) for line in result.stdout.splitlines())
            assert int(fields['maximum']) == literal
            assert int(fields['classes']) == 3**N
            witness = sum(int(x) << (2*k) for k,x in enumerate(fields['witness']))
            assert score(witness) == literal
            outputs.append(fields['witness'])
        assert outputs[0] == outputs[1]
        reports.append({'N':N,'m':m,'depth':depth,'flip_cap':cap,'free_guard':free,
                        'shells':len(chosen),'literal_centres':4**N,'maximum':str(literal),
                        'reference_and_grouped_agree':True})
    print(json.dumps({'precomputation_calibration':'PASS','cases':reports}),flush=True)
    return reports
