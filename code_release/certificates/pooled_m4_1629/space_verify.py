import sys

# The distance routine below uses int.bit_count(), which arrives in Python 3.10. macOS still ships
# 3.9 as `python3`, and the appendix tells a referee to run `python3 -I space_verify.py <index>`, so
# without this guard the reproducibility section's own command ends in an AttributeError traceback on
# a stock Mac. A reader who sees a stack trace has no reason to believe the certificate either, which
# is the whole point of shipping it. Refuse with an instruction instead.
if sys.version_info < (3, 10):
    sys.stderr.write(
        "NOT RUN: this verifier needs Python 3.10 or newer; it uses int.bit_count().\n"
        f"This interpreter is {sys.version.split()[0]}, and on macOS `python3` is usually 3.9.\n"
        "Run it with a newer interpreter, for example: python3.11 -I space_verify.py <index>\n")
    raise SystemExit(2)

import base64, hashlib, itertools, json, resource, time
from fractions import Fraction

started = time.perf_counter()
with open(sys.argv[1], "rb") as handle:
    source = handle.read()
index = json.loads(source)
assert all(len(k) == 22 and k.endswith("NNN") for k in index)
P = sorted(k[:19] for k in index)
assert len(P) == len(set(P)) == 1629
digest = lambda rows: hashlib.sha256(("\n".join(rows) + "\n").encode()).hexdigest()
assert digest(P) == "b749f0869f5d4a97d1e6c6c48f750694dbec54460846ad4730b5cba76136fa02"

payload = base64.b64decode("""
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 i4i1BTfe
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
TwAocfcQNlPhvIpDEE dc+aaAJHHhZtWbdq5H1CM4qefgNNJRmEcXQ+lQQZ45urIjcjTzREU8PSAUlo1XJiRwsyI1XMBIDkS/aQdKU2sg4d0RXcTjJjeiN2YP
VKiJrAIVUC5M5GgKuHFMwm+Igz1D8/sBjYmLvkLljfD0kSayfZMssw64ESOpaLwhkhfopv4CCnZJPzqWQ3VjorocLBqaFEVY9CL4SOnAlxeDFO2CNGRJEO5G
C8gMc4Y1if2vFynE+9nR4Levc3MXG5uQKqia3QnQhQfAOuO/ZVUkbSxRR6ULB5Q/SE3+gKxzjxG4WySHWmGrWHPUGcWorUCZFMiikxImrQCXBvgFDiTQQ7MB
sn+lF3jb2NxXQMZckS13GdLuXMVDRRJWKBUewx9uioPVZ9/iXkpGiB34N9/IBByWsfoUZCnyQPjlLtibv+X9PC+zx1ErIghUK+8QdfYAEA1U8xBZB8RiSYpA
YUZPBvyTioakglJejEsEowKjDFfIIOeYzoamXnguKjPkg6yRhGOYa2Q2cj2KjrVU96BAWg8FQ/RmDEPuMI4v+aM9lDh2ucND3VuaJ0aDwDFlR/5DxhHM7kMI
Q46pgSt3C2V6UAbrTDbY+h44LaQRS9Frg2Khuf2Bz8lJMcLQf+qgDVKnE2lYIrqIvDGE5bBwcE1hwFWq/88MruyQ6ITKjcuEZao8PwmlY6GHO/GP3lkNZ0eX
xT2WPi7UEWiAEpdpaZhkapO/UbWdIxsFntxEEsaFPemkWb5WwbPnCDMn2IRKJEj5ZMshMIyDSjO7QR0aaXWiXbIYCSU7ieR5A/iGB1M4ZCr1kPGg0R5+7CFZ
jmQlkrvVxQwA9+W1OiPKBaERHHEdIENVOYTu5TCmRrgiq3uCx104tI02FyLBXLagxkcZLmEOSBJO3PZEC90UWWYd0kH+Ws73SjViOvz6FAkmUnjwPo8C1SSa
1BrXUD3ZMa6tmTHrvEdNIQl5w5lBXuFv9CgvMhJjevwekAj/oFIFPh0sqcz04L9Skmi6lgWf8UsMXSYPgyL1NMCGL/lK9ByrGvyzA5yQM89yJXg1qJMI9saw
Vs2O29R2Yghj3EySCouYlBdiuTnRzdKm6XpoIIvaQGKjESCNDLnolUBggFU0rAMsiip5iXq4314jWcDVatl1kQFOyqXvxL/rMJUDIdZVswQ/dYRL4IbBY4ZF
z7NLXJypE+qQX7p8NMVojUM9in+ECeqy0sn3M7OhECou+7wIdMiwdqSsKtCBXMGUsVL74EC5BzJMGUVu8G8Q3OmEcv6Y/OLZ+efHQ6mGuPwyhCJGu8C/i3mE
4mZO6z/MJd6YOCoETzm29gZlh1h8W8slPiKTQhvTCqLDw/QIdy6oiXwLGAKjtlwH0oiSEpdNgSkhCwhe+ZQyWL8MwNJXWqScq4MbrYh8+Hy5CnXT1jyG0qiC
7f2VrAppjR3wm9b2XmaPrAJVlL8Ab9sgBHZhvCCiECxqxdiNX3OQEAcP5wRRnpqNPBI+3wPDlF+A5YgFakbXAQwyB9iSFzy5TF6DcPx6UxiG231vXI1d58au
1DjdchnGHk1OzNy6W3s2OUdWErNUn4XxOVg9atA2xg7gJteXUcrbERnKGULhZPuLgMmoWb0x9Ziy7+rjyqxPKDj00L/mR3cij3d+EoxLYrnJs3kkrS6N6A==
""")
assert len(payload) == 7738
assert hashlib.sha256(payload).hexdigest() == "ea65262c382b49cd431a99ac909213366b644fa7cefe8c317890f2e08ecbe895"
number = int.from_bytes(payload, "big")
assert number & 3 == 0
combinations = list(itertools.combinations(range(19), 4))
W = []
for i in range(3258):
    code = (number >> (len(payload)*8 - 19*(i+1))) & ((1 << 19)-1)
    rank, ternary = divmod(code, 81)
    assert 0 <= rank < len(combinations)
    choices = [0]*4
    for j in range(3, -1, -1):
        ternary, choices[j] = divmod(ternary, 3)
    word = list(P[i//2])
    for position, choice in zip(combinations[rank], choices):
        word[position] = [c for c in "ACGT" if c != P[i//2][position]][choice]
    W.append("".join(word))
assert len(set(W)) == 3258
assert digest(W) == "7213bb411954d2247e6df1e0f31533f6bf88fa87499048e9e7f14e36fe052a6e"

encode = lambda word: sum("ACGT".index(c) << (2*k) for k, c in enumerate(word))
mask = sum(1 << (2*k) for k in range(19))
def distance(a, b):
    difference = a ^ b
    return ((difference | (difference >> 1)) & mask).bit_count()

EP, EW = list(map(encode, P)), list(map(encode, W))
for i, word in enumerate(EW):
    for j, pattern in enumerate(EP):
        d = distance(word, pattern)
        assert d == 4 if j == i//2 else d >= 5

higher = [set() for _ in W]
for i in range(len(W)):
    for j in range(i+1, len(W)):
        if distance(EW[i], EW[j]) <= 8:
            higher[i].add(j)

def common_ball(a, b, c):
    states = {(0, 0, 0)}
    for ca, cb, cc in zip(a, b, c):
        next_states = set()
        for da, db, dc in states:
            for letter in "ACGT":
                state = (da+(letter != ca), db+(letter != cb), dc+(letter != cc))
                if max(state) <= 4:
                    next_states.add(state)
        states = next_states
        if not states:
            return False
    return bool(states)

triples = 0
for i in range(len(W)):
    for j in higher[i]:
        for k in higher[i] & higher[j]:
            triples += 1
            assert not common_ball(W[i], W[j], W[k]), (i, j, k)
assert triples == 7404
assert common_ball(W[0], W[1], P[0])

print(json.dumps({
    "status": "PASS",
    "support": len(W),
    "mass": str(Fraction(len(W), 2)),
    "pair_checks": len(W)*(len(W)-1)//2,
    "triangles_checked": triples,
    "worst_ball_mass_bound": "1",
    "negative_control_detected": True,
    "source_sha256": hashlib.sha256(source).hexdigest(),
    "support_sha256": digest(W),
    "seconds": time.perf_counter()-started,
    "peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss /
                    (2**20 if sys.platform == "darwin" else 1024)
}))
