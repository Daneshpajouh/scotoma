"""Loading a guide panel.

A panel is a set of distinct guides, held as packed window codes. Only the spacer bits are ever
compared, so a panel may be given as 20-mer spacers or as full 23-mer windows and the two are
interchangeable.

Duplicate spacers are collapsed. That matters for a density: counting one guide twice would
double its contribution to the minimum-over-panel distance for no biological reason.
"""
from __future__ import annotations

import csv
import sys

import numpy as np

from .coding import K, SPACER_LEN, encode_spacer, encode_window, spacer_of


def _as_code(token: str):
    token = token.strip().upper()
    if len(token) == K:
        return encode_window(token)
    if len(token) == SPACER_LEN:
        return encode_spacer(token)
    return None


def load_panel_text(path):
    """Read a panel from a text file of one guide per line.

    Blank lines and lines beginning with ``#`` are ignored. Each remaining line must be a
    20-mer spacer or a 23-mer window of ACGT.

    Returns
    -------
    tuple
        ``(codes, report)`` where ``codes`` is a sorted uint64 array of distinct guides and
        ``report`` records how many lines were read, skipped and deduplicated.
    """
    seen = {}
    n_lines = 0
    rejected = []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            n_lines += 1
            code = _as_code(line)
            if code is None:
                rejected.append(line[:40])
                continue
            seen.setdefault(spacer_of(code), code)
    return _finish(seen, n_lines, rejected, source=str(path))


def load_panel_changeseq(path, spacer_column="sgRNA", window_column="Align.off-target",
                         mismatch_column="Align.#Mismatches", bulge_column="Align.#Bulges"):
    """Read a panel from a CHANGE-seq style measurement table.

    The guide set is taken from the rows that are each guide's own on-target, meaning zero
    mismatches and zero bulges. Those rows carry the genomic window the guide was designed
    against, so the panel comes from the same table as any measurements, with no separate
    design file to drift out of step.

    Rows whose mismatch or bulge field is missing or non-numeric are skipped rather than
    guessed at.
    """
    # Measurement tables carry very long alignment fields; raise the limit rather than crash.
    try:
        csv.field_size_limit(1 << 30)
    except OverflowError:  # pragma: no cover - platform dependent ceiling
        csv.field_size_limit(sys.maxsize)

    seen = {}
    n_rows = 0
    rejected = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                mm = int(row[mismatch_column] or 99)
                bulges = int(row[bulge_column] or 0)
            except (TypeError, ValueError, KeyError):
                continue
            if mm != 0 or bulges != 0:
                continue
            n_rows += 1
            window = (row.get(window_column) or "").upper().replace("-", "")
            spacer = (row.get(spacer_column) or "").upper()[:SPACER_LEN]
            code = encode_window(window) if len(window) == K else None
            if code is None and len(spacer) == SPACER_LEN:
                code = encode_spacer(spacer)
            if code is None:
                rejected.append((spacer or window)[:40])
                continue
            seen.setdefault(spacer_of(code), code)
    return _finish(seen, n_rows, rejected, source=str(path))


def _finish(seen, n_read, rejected, source):
    codes = np.array(sorted(seen.values()), dtype=np.uint64)
    report = {
        "source": source,
        "rows_considered": n_read,
        "distinct_guides": int(codes.size),
        "duplicates_collapsed": max(0, n_read - len(rejected) - int(codes.size)),
        "rejected": len(rejected),
        "rejected_examples": rejected[:5],
    }
    return codes, report


def load_panel(path, changeseq: bool = False):
    """Dispatch to the text or CHANGE-seq reader.

    ``changeseq`` is explicit rather than sniffed from the extension, because silently
    misreading a panel would change every number downstream.
    """
    return load_panel_changeseq(path) if changeseq else load_panel_text(path)


def panel_from_windows(windows, n: int, seed: int = 0, gc_match_to=None,
                       unique_only: bool = True, max_multiplicity: int = 2):
    """Draw a synthetic panel of ``n`` guides from a pool of oriented window codes.

    This is the discriminating control for a matched-arm density comparison. A
    composition-preserving shuffle (see :func:`scotoma.coding.shuffle_panel`) tests whether a result
    needs the *real* guide sequences, but a shuffled spacer is not a genuine genomic 20-mer and so
    cannot separate two different explanations: that a clinical panel's near-cognate sites are
    reference-present because the panel was developed against that reference, against the weaker
    claim that any real genomic guide resembles typical genomic sequence more than atypical
    sequence. Drawing a panel from real reference-present windows holds "genuine genomic sequence"
    fixed and varies only which guides.

    ``gc_match_to`` optionally supplies an existing panel whose spacer-GC histogram the drawn panel
    must match bin by bin, so the drawn panel cannot differ from the clinical one in composition.
    Clinical guides are substantially GC-richer than random genomic protospacers, so without this
    the control would be confounded by exactly the offset that invalidated four earlier null models.
    """
    import numpy as _np
    from .coding import spacer_gc

    windows = _np.asarray(windows, dtype=_np.uint64)
    rng = _np.random.default_rng(seed + 31000)

    # Draw from DISTINCT codes, not from window instances. Drawing uniformly over instances
    # samples each code in proportion to how often it occurs, so a code present ten thousand
    # times is ten thousand times likelier to be picked. A first version of this control did
    # exactly that and produced a panel of highly repetitive sequence: 110 guides whose codes
    # together accounted for 782,752 window instances, about 7,100 copies each, where a clinical
    # panel's guides are essentially single-site. GC was matched and multiplicity was not, so the
    # control answered a question nobody asked.
    multiplicity_report = {}
    if unique_only:
        uniq, counts = _np.unique(windows, return_counts=True)
        multiplicity_report = {
            "distinct_codes_in_pool": int(uniq.size),
            "max_multiplicity_in_pool": int(counts.max()) if counts.size else 0,
            "max_multiplicity_allowed": int(max_multiplicity),
        }
        eligible = uniq[counts <= max_multiplicity]
        multiplicity_report["eligible_codes"] = int(eligible.size)
        if eligible.size < n:
            raise ValueError(
                f"only {eligible.size} codes have multiplicity <= {max_multiplicity}; "
                f"{n} guides requested"
            )
        windows = eligible

    if windows.size < n:
        raise ValueError(f"pool of {windows.size} is smaller than the {n} guides requested")

    if gc_match_to is None:
        idx = rng.choice(windows.size, n, replace=False)
        drawn = windows[idx]
        report = {"n": int(drawn.size), "gc_matched": False, "pool": int(windows.size)}
        report.update(multiplicity_report)
        return drawn, report

    target_gc = spacer_gc(_np.asarray(gc_match_to, dtype=_np.uint64))
    pool_gc = spacer_gc(windows)
    picks = []
    shortfall = {}
    for b in range(int(target_gc.max()) + 1):
        want = int((target_gc == b).sum())
        if want == 0:
            continue
        avail = _np.flatnonzero(pool_gc == b)
        take = min(want, avail.size)
        if take < want:
            shortfall[b] = want - take
        if take:
            picks.append(rng.choice(avail, take, replace=False))
    if not picks:
        raise ValueError("no pool window matched any GC bin of the reference panel")
    drawn = windows[_np.concatenate(picks)]
    drawn_gc = spacer_gc(drawn)
    max_diff = max(abs(int((target_gc == b).sum()) - int((drawn_gc == b).sum()))
                   for b in range(21))
    report = {
        "n": int(drawn.size),
        "gc_matched": True,
        "pool": int(windows.size),
        "requested": int(target_gc.size),
        "gc_shortfall_by_bin": shortfall,
        "max_per_bin_difference": int(max_diff),
        "mean_spacer_gc_target": float(target_gc.mean()),
        "mean_spacer_gc_drawn": float(drawn_gc.mean()),
    }
    report.update(multiplicity_report)
    return drawn, report
