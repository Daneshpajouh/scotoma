"""Exhaust the ACGTA example: unrestricted K=4, minimum sound cover=5."""
from itertools import product
import json


def main():
    guide = (0,1,2,3,0)
    library = sorted({guide[:i]+guide[i+1:] for i in range(len(guide))})
    words = list(product(range(4), repeat=4))
    balls = [sum(1 << i for i,x in enumerate(words) if sum(a != b for a,b in zip(x,c)) <= 2)
             for c in words]
    universe = 0
    for p in library:
        universe |= balls[words.index(p)]
    covered = [b & universe for b in balls]
    sound = [i for i,b in enumerate(balls) if b & ~universe == 0]

    def search(indices, limit):
        failed = set()
        owners = {i:[j for j in indices if covered[j] & (1 << i)]
                  for i in range(len(words)) if universe & (1 << i)}
        def recurse(missing, slots):
            if missing == 0:
                return []
            if slots == 0 or (missing,slots) in failed:
                return None
            best = max((covered[j] & missing).bit_count() for j in indices)
            if best*slots < missing.bit_count():
                return None
            point = min((i for i in owners if missing & (1 << i)), key=lambda i:len(owners[i]))
            candidates = sorted(owners[point], key=lambda j:-(covered[j] & missing).bit_count())
            for j in candidates:
                tail = recurse(missing & ~covered[j],slots-1)
                if tail is not None:
                    return [j]+tail
            failed.add((missing,slots))
            return None
        return recurse(universe,limit),len(failed)

    no_three, states = search(list(range(len(words))),3)
    four, _ = search(list(range(len(words))),4)
    no_sound_four, sound_states = search(sound,4)
    sound_five, _ = search(sound,5)
    assert len(library) == 5 and no_three is None and four is not None and len(four) == 4
    assert no_sound_four is None and sound_five is not None and len(sound_five) == 5
    letters='ACGT'
    show=lambda i:''.join(letters[x] for x in words[i])
    print(json.dumps({'passed':True,'guide':'ACGTA','radius':2,'ambient_words':len(words),
          'union_words':universe.bit_count(),'library_size':len(library),'K':4,'minimum_sound_cover':5,
          'four_cover':[show(i) for i in four],'sound_cover':[show(i) for i in sound_five],
          'sound_centres':len(sound),'failed_three_slot_states':states,
          'failed_sound_four_slot_states':sound_states}),flush=True)


if __name__ == '__main__':
    main()
