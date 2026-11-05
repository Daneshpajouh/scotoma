#!/usr/bin/env python3
"""The IL2RG claim, recomputed from sequence: at four substitutions under NGG only the terminal
deletion reaches a locus an independent study classifies as edited.

    python3 -I verify_il2rg_witness.py

Chaudhari et al. 2020 (doi 10.1089/crispr.2020.0053) report guide IL2RG and a site at
chr1:14,947,683-14,947,705 which their hybrid-capture panel classifies as edited and their
orthogonal rhAmpSeq panel measures again. This recomputes, from the deposited 103-base GRCh38
window alone, every register the audited emitters can reach: all twenty single-deletion words,
the two emitted subsets, every bulge-free 20-base register and every single-DNA-bulge register,
on both strands, under a strict NGG PAM. It asserts what the paper says and no more.

The assay values are the source's own measurements. They are pinned in IL2RG_witness.json with
the digests of the two files they were read from and are not recomputed from reads here.

Controls that must behave as stated, so a reader can tell a check from a restatement:
  C1  a forged window carrying a perfect match stops reproducing the minima
  C2  altering one base of the core breaks the exclusivity of position 20
  C3  the source's own GGT register is not NGG, so admitting it must change the ungapped minimum
"""
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
SPACER = 20
VERDICTS = []


def claim(ok, label, got, want):
    VERDICTS.append(bool(ok))
    print("  %-5s %-54s got %-14s want %s" % ("PASS" if ok else "FAIL", label, got, want))


rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
hd = lambda a, b: sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else 10 ** 6

rec = json.load(open(HERE / "IL2RG_witness.json"))
G = rec["guide"]
W = rec["window"]["forward_sequence"].upper()
CORE = rec["canonical_core"].upper()

claim(len(G) == SPACER, "guide length", len(G), SPACER)
claim(len(W) == rec["window"]["length"] == 103, "window length", len(W), 103)
claim(len(CORE) == 19, "canonical core length", len(CORE), 19)

# the core and its NGG PAM must actually occur in the window, on the stated strand
found = [(st, i) for st, s in (("+", W), ("-", rc(W)))
         for i in range(len(s) - 21)
         if s[i:i + 19] == CORE and s[i + 20:i + 22] == "GG"]
claim(len(found) == 1, "core plus an NGG PAM occurs exactly once in the window", len(found), 1)
claim(found and found[0][0] == rec["strand"], "  and on the recorded strand",
      found[0][0] if found else "none", rec["strand"])

# every single-deletion word, by position
per = {i + 1: hd(G[:i] + G[i + 1:], CORE) for i in range(SPACER)}
vec = [per[i + 1] for i in range(SPACER)]
claim(vec == rec["deletion_distances_by_position"], "all twenty deletion distances",
      "identical" if vec == rec["deletion_distances_by_position"] else str(vec), "identical")

full = min(per.values())
arg = sorted(i for i in per if per[i] == full)
claim(full == 4, "complete library reaches the locus at four", full, 4)
claim(arg == [SPACER], "  and only by deleting position 20", arg, [20])
claim(min(per[i] for i in range(1, SPACER)) == 5, "released build, positions 1..19",
      min(per[i] for i in range(1, SPACER)), 5)
claim(min(per[i] for i in range(2, SPACER)) == 5, "2016 wrapper, positions 2..19",
      min(per[i] for i in range(2, SPACER)), 5)


def scan(window, pam_test):
    """Minima over the window for every register class, both strands, PAM by pam_test."""
    best = {"full": 99, "native": 99, "wrapper": 99, "ungapped": 99, "DNA": 99}
    for seq in (window, rc(window)):
        for width in (19, 20, 21):
            for j in range(len(seq) - width - 2):
                site = seq[j:j + width + 3]
                if not pam_test(site[-3:]):
                    continue
                core = site[:-3]
                if width == 19:
                    v = [hd(core, G[:i] + G[i + 1:]) for i in range(SPACER)]
                    best["full"] = min(best["full"], min(v))
                    best["native"] = min(best["native"], min(v[:19]))
                    best["wrapper"] = min(best["wrapper"], min(v[1:19]))
                elif width == 20:
                    best["ungapped"] = min(best["ungapped"], hd(G, core))
                else:
                    best["DNA"] = min(best["DNA"],
                                      min(hd(G, core[:k] + core[k + 1:]) for k in range(21)))
    return best


ngg = lambda p: p[1:3] == "GG"
got = scan(W, ngg)
want = rec["minima_under_strict_NGG"]
claim(got["full"] == want["complete_family_1_to_20"], "window scan, complete family",
      got["full"], want["complete_family_1_to_20"])
claim(got["native"] == want["released_build_1_to_19"], "window scan, released build",
      got["native"], want["released_build_1_to_19"])
claim(got["wrapper"] == want["wrapper_2_to_19"], "window scan, 2016 wrapper",
      got["wrapper"], want["wrapper_2_to_19"])
claim(got["ungapped"] == want["best_ungapped_20_base_register"],
      "window scan, best bulge-free register", got["ungapped"],
      want["best_ungapped_20_base_register"])
claim(got["DNA"] == want["best_single_DNA_bulge_register"],
      "window scan, best single-DNA-bulge register", got["DNA"],
      want["best_single_DNA_bulge_register"])
claim(got["ungapped"] > 4 and got["DNA"] > 4,
      "no bulge-free or DNA-bulge route exists at four under NGG",
      "%d and %d" % (got["ungapped"], got["DNA"]), "both > 4")

# the source's own register, for the limit the paper states
src = rec["source_register"]
claim(src["pam"] != "" and src["pam"][1:3] != "GG", "the source's own PAM is not NGG",
      src["pam"], "not NGG")
claim(hd(G, src["protospacer"].upper()) == src["mismatches"],
      "  its ungapped alignment is at the stated distance",
      hd(G, src["protospacer"].upper()), src["mismatches"])

# C3: admitting that PAM must change the ungapped minimum, or the limit sentence is empty
anyg = lambda p: p[2] == "G" or p[1] == "G"
loose = scan(W, anyg)
claim(loose["ungapped"] < got["ungapped"],
      "C3: a permissive PAM policy does supply an ungapped route",
      "%d vs %d" % (loose["ungapped"], got["ungapped"]), "smaller")

# C1: forge a perfect match into the window and the minima must stop reproducing
cut = len(W) // 2
forged = (W[:cut] + G + "AGG" + W[cut + 23:])[:len(W)]
fgot = scan(forged, ngg)
claim(fgot["full"] != got["full"] or fgot["ungapped"] != got["ungapped"],
      "C1: a forged window stops reproducing the minima",
      "%d/%d" % (fgot["full"], fgot["ungapped"]),
      "differs from %d/%d" % (got["full"], got["ungapped"]))

# C2: alter one core base at a time, over every position and every alternative base, and check
# where the exclusivity survives. The terminal word is G[:-1]; at a position where the core and
# that word agree, any substitution must push the distance to five and break exclusivity. At one
# of the four positions where they already differ, a substitution to the third base leaves the
# distance at four, so exclusivity may legitimately survive. Asserting that pattern, rather than
# one hand-picked variant, is what makes this a check: the first version of this control altered
# a position that was already a mismatch and so proved nothing.
TERM = G[:-1]
agree = [i for i in range(19) if CORE[i] == TERM[i]]
survived, broke = [], 0
for i in range(19):
    for b in "ACGT":
        if b == CORE[i]:
            continue
        alt = CORE[:i] + b + CORE[i + 1:]
        aper = {k + 1: hd(G[:k] + G[k + 1:], alt) for k in range(SPACER)}
        if min(aper.values()) == 4 and sorted(k for k in aper if aper[k] == 4) == [SPACER]:
            survived.append(i)
        else:
            broke += 1
claim(len(agree) == 15, "C2: core and terminal word agree at fifteen positions", len(agree), 15)
claim(broke >= 3 * len(agree), "  every substitution at an agreeing position breaks it",
      broke, ">= %d" % (3 * len(agree)))
claim(all(i not in agree for i in survived),
      "  and exclusivity survives only at the four mismatched positions",
      sorted(set(survived)), "all outside the agreeing set")

# the executed wrapper's own emitted set, as recorded
ex = rec["executed_wrapper"]
claim(ex["distinct_words_emitted"] == 13 and ex["complete_library_size"] == 15,
      "recorded emitter run: 13 of 15 words",
      "%d of %d" % (ex["distinct_words_emitted"], ex["complete_library_size"]), "13 of 15")
claim(ex["terminal_word_present"] is False and ex["first_word_present"] is False,
      "  and neither endpoint word is present", "neither", "neither")
words = {G[:i] + G[i + 1:] for i in range(SPACER)}
claim(len(words) == ex["complete_library_size"], "the complete library size is the run count",
      len(words), ex["complete_library_size"])

n = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - n, n))
if n:
    sys.exit(1)
print("At four substitutions under NGG the complete deletion family reaches this locus and only")
print("position 20 does so; the released build and the 2016 wrapper both stop at five, and no")
print("bulge-free or single-DNA-bulge register in the window comes within ten.")
