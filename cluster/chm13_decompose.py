#!/usr/bin/env python3
"""Decompose the T2T-CHM13 candidate gain into newly represented sequence and copy-number change.

The reviewer's objection to reporting a bare "assembly effect" is that it is a composite: a guide
can gain candidates because CHM13 contains sequence GRCh38 omits, or because sequence present in
both is represented at a different copy number. Those have different meanings and the paper should
not merge them.

Both components are computable from the two universes alone, with no cross-assembly alignment. For
each guide, take its CHM13 candidates at a threshold and ask whether each candidate's canonical code
is also in the GRCh38 universe:

  shared  : code present in GRCh38 too, so any count difference is copy-number representation
  novel   : code absent from GRCh38, so the candidate lives in sequence GRCh38 does not contain

Counts are multiplicity-summed genomic sites under the same convention as the main results.
"""
import argparse, json, time

import numpy as np

from scotoma.coding import revcomp64
from scotoma.distance import spacer_mismatches
from scotoma.panel import load_panel
from scotoma.universe import load_universe

GG = np.uint64(0b1010)
LAST2 = np.uint64(0b1111)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chm13", required=True)
    ap.add_argument("--grch38", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--changeseq", action="store_true")
    ap.add_argument("--thresholds", type=int, nargs="+", default=[3, 4, 5, 6])
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    t0 = time.time()
    c_codes, c_mult, c_meta = load_universe(a.chm13, with_mult=True)
    print(f"CHM13 {c_codes.size:,} codes, {time.time()-t0:.0f}s", flush=True)
    g_codes, _ = load_universe(a.grch38)
    print(f"GRCh38 {g_codes.size:,} codes, {time.time()-t0:.0f}s", flush=True)
    guides, rep = load_panel(a.panel, changeseq=a.changeseq)
    print(f"panel {rep['distinct_guides']} guides", flush=True)

    ths = tuple(a.thresholds)
    shared = [{t: 0 for t in ths} for _ in range(guides.size)]
    novel = [{t: 0 for t in ths} for _ in range(guides.size)]

    n_chunks = (c_codes.size + a.chunk - 1) // a.chunk
    for ci, start in enumerate(range(0, c_codes.size, a.chunk)):
        stop = min(start + a.chunk, c_codes.size)
        block = c_codes[start:stop]
        rc = revcomp64(block)
        ok_f = (block & LAST2) == GG
        ok_r = (rc & LAST2) == GG
        mult = c_mult[start:stop].astype(np.int64)
        # membership of this block's codes in GRCh38, computed once per chunk
        idx = np.clip(np.searchsorted(g_codes, block), 0, g_codes.size - 1)
        in_g38 = g_codes[idx] == block
        far = np.int8(99)
        for gi in range(guides.size):
            g = int(guides[gi])
            best = np.where(ok_f, spacer_mismatches(block, g), far)
            np.minimum(best, np.where(ok_r, spacer_mismatches(rc, g), far), out=best)
            for t in ths:
                hit = best <= t
                if not hit.any():
                    continue
                shared[gi][t] += int(mult[hit & in_g38].sum())
                novel[gi][t] += int(mult[hit & ~in_g38].sum())
        if (ci + 1) % 10 == 0 or ci + 1 == n_chunks:
            print(f"  chunk {ci+1}/{n_chunks}  {time.time()-t0:.0f}s", flush=True)

    rows = [{"guide_index": i, "shared": shared[i], "novel": novel[i]} for i in range(guides.size)]
    payload = {"chm13": c_meta, "thresholds": list(ths), "guides": int(guides.size), "rows": rows,
               "wall_seconds": round(time.time() - t0, 1)}
    with open(a.out_json, "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"wrote {a.out_json}", flush=True)

    for t in ths:
        sh = [r["shared"][t] for r in rows]
        nv = [r["novel"][t] for r in rows]
        tot = [s + n for s, n in zip(sh, nv)]
        frac = sorted(n / x for n, x in zip(nv, tot) if x)
        med = frac[len(frac) // 2]
        print(f"<= {t} mm: novel-sequence share of CHM13 candidates, "
              f"min {min(frac):.4f} median {med:.4f} max {max(frac):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
