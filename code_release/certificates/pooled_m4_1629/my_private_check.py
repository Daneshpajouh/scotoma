# Independent: library rebuilt from the guide list (not the census file); plain distance; no shared code.
#
# Three things were wrong for a reader and are fixed here, none of them touching the certificate.
# The guide list was a required argument that ships only as JSON, so the appendix's own command,
# `python3 -I my_private_check.py` with no argument, raised IndexError instead of printing the
# transcript beneath it; the panel now loads by default from the released panel reference. The two
# certificate files were opened by bare name, so the script ran only from inside this directory;
# they are now resolved beside the script. And a missing input exited 1 through a traceback, which
# reads as a refuted claim rather than an absent file; it now exits 2 and says what to supply.
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
# The RELEASED layout first, the working-tree sibling second, the same order verify_minimality.py
# uses. This file is copied into two trees whose depths differ, so a single relative path is correct
# in one of them and silently wrong in the other.
PANEL_CANDIDATES = [
    os.path.join(HERE, "..", "..", "data", "producers_2026-10-20", "PANELREF_GRCh38_NGG.json"),
    os.path.join(HERE, "..", "..", "producers_2026-10-20", "PANELREF_GRCh38_NGG.json"),
]
PANEL = next((p for p in PANEL_CANDIDATES if os.path.exists(p)), PANEL_CANDIDATES[0])


def fail_missing(path, what):
    sys.stderr.write(f"NOT RUN: {what} not found at {path}\n"
                     "This is a missing input, not a failed check. Supply the file, or pass a\n"
                     "guide list as the first argument: one guide per line, guide in the last column.\n")
    raise SystemExit(2)


source = sys.argv[1] if len(sys.argv) > 1 else PANEL
if not os.path.exists(source):
    fail_missing(source, "the guide list")
if source.endswith(".json"):
    # The released panel reference carries the 110 spacers as rows[i]["spacer"].
    guides = [r["spacer"] for r in json.load(open(source))["rows"]]
else:
    guides = [l.split()[-1].strip() for l in open(source) if l.strip()]

lib = sorted({g[:j] + g[j+1:] for g in guides for j in range(20)})
for name in ("P_1629.txt", "W_3258.txt"):
    if not os.path.exists(os.path.join(HERE, name)):
        fail_missing(os.path.join(HERE, name), name)
P = [l.strip() for l in open(os.path.join(HERE, "P_1629.txt")) if l.strip()]
W = [l.strip() for l in open(os.path.join(HERE, "W_3258.txt")) if l.strip()]
same = lib == sorted(P)
print("library rebuilt from guides:", len(lib), "| equals census library:", same)
ham = lambda a, b: sum(x != y for x, y in zip(a, b))
bad_in_U = bad_private = 0
owners = {}
for w in W:
    ds = [ham(w, q) for q in lib]
    close = [i for i, d in enumerate(ds) if d <= 4]
    if not close: bad_in_U += 1
    if len(close) != 1: bad_private += 1
    else: owners[close[0]] = owners.get(close[0], 0) + 1
two = sum(1 for v in owners.values() if v == 2)
print("words:", len(W), "distinct:", len(set(W)), "| outside U:", bad_in_U, "| not private:", bad_private,
      "| library words owning exactly two:", two, "of", len(lib))
# Exit status so the run can gate: every word private and inside U, every pattern owning exactly two.
ok = same and bad_in_U == 0 and bad_private == 0 and two == len(lib) == len(P)
print("PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
