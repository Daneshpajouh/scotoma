"""The reference universe: which SpCas9 target windows a reference genome contains.

A universe is a sorted array of distinct **canonical** codes plus the multiplicity of each. Two
things follow from storing it this way. Membership becomes a binary search rather than a rescan,
so partitioning a haplotype against a reference costs one pass over the haplotype. And
multiplicity is a per-guide site count across both strands, because canonicalisation folds the
two strands together.

Build one per reference you want to test against, typically GRCh38 and T2T-CHM13. The file
format is a plain ``.npz`` with ``codes``, ``mult`` and a one-element ``meta`` array holding
JSON, so anything that can read numpy can read it without this package.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from .coding import DEFAULT_PAM, K, PAM_SCOPES, canonical, iter_fasta, oriented_sites


def build_universe(fasta_paths, label: str, primary_only: bool = False,
                   pam: str = DEFAULT_PAM):
    """Scan one or more FASTA files and return ``(codes, mult, meta)``.

    Parameters
    ----------
    fasta_paths : iterable of str
        FASTA files, gzipped or not. All records in all files are pooled into one universe.
    label : str
        Name recorded in the metadata, for example ``"GRCh38"``.
    primary_only : bool
        Skip records whose name contains ``_`` or ``.``, which excludes GRCh38 alts, patches
        and unplaced scaffolds under the usual naming. Off by default because it is a
        convention, not a guarantee, and silently dropping sequence from a reference universe
        would inflate every downstream absence count.

    Returns
    -------
    tuple
        ``codes`` sorted unique canonical uint64, ``mult`` int32 of the same length, and a
        ``meta`` dict.
    """
    if pam not in PAM_SCOPES:
        raise ValueError(f"pam must be one of {PAM_SCOPES}, got {pam!r}")
    started = time.time()
    parts = []
    total_sites = 0
    n_records = 0
    assembled_bp = 0
    skipped = []
    for path in fasta_paths:
        for name, bases in iter_fasta(path):
            if primary_only and ("_" in name or "." in name):
                skipped.append(name)
                continue
            if bases.size < K:
                continue
            n_records += 1
            assembled_bp += int(bases.size)
            sites = oriented_sites(bases, pam=pam)
            if sites.size:
                parts.append(canonical(sites))
                total_sites += int(sites.size)
    if not parts:
        raise ValueError("no usable sequence found; check the FASTA paths")

    codes = np.concatenate(parts)
    del parts
    codes.sort()
    edges = np.flatnonzero(np.concatenate(([True], codes[1:] != codes[:-1])))
    mult = np.diff(np.concatenate((edges, [codes.size]))).astype(np.int32)
    uniq = np.ascontiguousarray(codes[edges])

    meta = {
        "label": label,
        "k": K,
        "pam": pam,
        "records": n_records,
        "assembled_bp": assembled_bp,
        "target_sites_total": total_sites,
        "distinct_canonical": int(uniq.size),
        "single_site": int((mult == 1).sum()),
        "max_sites": int(mult.max()) if mult.size else 0,
        "primary_only": bool(primary_only),
        "skipped_records": len(skipped),
        "wall_seconds": round(time.time() - started, 1),
    }
    return uniq, mult, meta


def save_universe(path, codes, mult, meta) -> None:
    """Write a universe to ``.npz``, creating parent directories as needed."""
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    np.savez_compressed(
        path,
        codes=np.asarray(codes, dtype=np.uint64),
        mult=np.asarray(mult, dtype=np.int32),
        meta=np.array([json.dumps(meta)], dtype=object),
    )


def _read_meta(archive):
    """Recover a universe's metadata dict, whatever an older writer stored it as.

    Universes built by earlier tooling store the JSON as ``bytes`` rather than ``str``, and
    ``str(b'{...}')`` is ``"b'{...}'"``, which is not JSON: it fails at character 1. A real
    GS_CHM13 file hit exactly that. Metadata is descriptive and nothing in the computation depends
    on it, so an unreadable block degrades to an empty dict with a recorded reason instead of
    killing a job that had already loaded 226 million codes.
    """
    if "meta" not in archive:
        return {}
    try:
        raw = archive["meta"]
        item = raw[0] if getattr(raw, "shape", ()) else raw.item()
        if isinstance(item, bytes):
            item = item.decode("utf-8", errors="replace")
        elif not isinstance(item, str):
            item = str(item)
        value = json.loads(item)
        return value if isinstance(value, dict) else {"meta_raw": value}
    except Exception as exc:  # noqa: BLE001 - metadata must never break a load
        return {"meta_unreadable": f"{type(exc).__name__}: {exc}"}


def load_universe(path, with_mult: bool = False):
    """Load a universe written by :func:`save_universe`.

    Returns ``(codes, meta)``, or ``(codes, mult, meta)`` when ``with_mult`` is true.
    Multiplicities are a few hundred megabytes for a human reference, so they are left on disk
    unless asked for.

    Raises
    ------
    ValueError
        If the stored codes are not sorted and unique, which every membership test here
        assumes.
    """
    with np.load(path, allow_pickle=True) as z:
        codes = z["codes"]
        meta = _read_meta(z)
        mult = z["mult"] if (with_mult and "mult" in z) else None
    if codes.size > 1 and not (codes[1:] > codes[:-1]).all():
        raise ValueError(f"{path}: codes are not sorted and unique")
    if with_mult:
        return codes, mult, meta
    return codes, meta


def contains(universe_codes, query_codes):
    """Membership of each query in the universe, as a boolean mask.

    ``query_codes`` may be oriented; they are canonicalised here, which is the only place
    canonicalisation belongs. ``universe_codes`` must be sorted and unique, as
    :func:`load_universe` enforces.
    """
    universe_codes = np.asarray(universe_codes, dtype=np.uint64)
    q = canonical(query_codes)
    if universe_codes.size == 0:
        return np.zeros(q.size, dtype=bool)
    idx = np.searchsorted(universe_codes, q)
    idx_clipped = np.clip(idx, 0, universe_codes.size - 1)
    return universe_codes[idx_clipped] == q


def reference_candidate_counts(universe_codes, universe_mult, guide_code,
                               thresholds=(0, 1, 2, 3, 4, 5, 6), chunk: int = 4_000_000):
    """Candidates within each mismatch threshold of one guide, counted **in the reference**.

    This is what a conventional reference-anchored nominator reports, so it is the incumbent
    quantity: the `N` a validation campaign would be sized against.

    Orientation has to be recovered before distances are taken. A universe stores canonical codes,
    which have lost the strand the window was found on. Every stored window had a valid PAM in at
    least one orientation, so for each canonical code ``c`` the admissible oriented forms are
    whichever of ``c`` and ``revcomp(c)`` end in ``GG``, and the window's distance to the guide is
    the minimum over those. Taking the distance on the canonical code instead would silently
    impose a ``CC`` prefix on about half the reference.

    Returns
    -------
    dict
        ``{"distinct": {t: count}, "sites": {t: count}}``. ``distinct`` counts canonical windows;
        ``sites`` sums their both-strand multiplicities, which is the genomic-position count.
    """
    from .coding import K, revcomp64
    from .distance import spacer_mismatches

    universe_codes = np.asarray(universe_codes, dtype=np.uint64)
    universe_mult = (None if universe_mult is None
                     else np.asarray(universe_mult, dtype=np.int64))
    thresholds = tuple(int(t) for t in thresholds)
    distinct = {t: 0 for t in thresholds}
    sites = {t: 0 for t in thresholds}

    # Last two bases are bits 3..0; GG is 0b1010.
    gg = np.uint64(0b1010)
    last_two = np.uint64(0b1111)

    for start in range(0, universe_codes.size, chunk):
        stop = min(start + chunk, universe_codes.size)
        block = universe_codes[start:stop]
        rc = revcomp64(block)
        ok_f = (block & last_two) == gg
        ok_r = (rc & last_two) == gg

        best = np.full(block.size, 99, np.int8)
        df = spacer_mismatches(block, guide_code)
        np.minimum(best, np.where(ok_f, df, np.int8(99)), out=best)
        dr = spacer_mismatches(rc, guide_code)
        np.minimum(best, np.where(ok_r, dr, np.int8(99)), out=best)

        mult = (universe_mult[start:stop] if universe_mult is not None
                else np.ones(block.size, np.int64))
        for t in thresholds:
            hit = best <= t
            distinct[t] += int(hit.sum())
            sites[t] += int(mult[hit].sum())

    return {"distinct": distinct, "sites": sites}


def panel_reference_counts(universe_codes, universe_mult, guide_codes,
                           thresholds=(0, 1, 2, 3, 4, 5, 6), chunk: int = 2_000_000,
                           progress=None):
    """:func:`reference_candidate_counts` for a whole panel, in one pass over the universe.

    Chunking is the outer loop and the panel the inner one, so each chunk's reverse complements
    and PAM admissibility are computed once and reused across every guide. For a 110-guide panel
    that is 110 times less reverse-complement work than calling the single-guide function in a
    loop.

    Returns
    -------
    list of dict
        One entry per guide, in the order given, each ``{"guide_index", "distinct", "sites"}``.
    """
    from .coding import revcomp64
    from .distance import spacer_mismatches

    universe_codes = np.asarray(universe_codes, dtype=np.uint64)
    universe_mult = (None if universe_mult is None
                     else np.asarray(universe_mult, dtype=np.int64))
    guides = [int(g) for g in np.asarray(guide_codes, dtype=np.uint64).tolist()]
    thresholds = tuple(int(t) for t in thresholds)

    distinct = [{t: 0 for t in thresholds} for _ in guides]
    sites = [{t: 0 for t in thresholds} for _ in guides]

    gg = np.uint64(0b1010)
    last_two = np.uint64(0b1111)
    n_chunks = (universe_codes.size + chunk - 1) // chunk

    for ci, start in enumerate(range(0, universe_codes.size, chunk)):
        stop = min(start + chunk, universe_codes.size)
        block = universe_codes[start:stop]
        rc = revcomp64(block)
        ok_f = (block & last_two) == gg
        ok_r = (rc & last_two) == gg
        mult = (universe_mult[start:stop] if universe_mult is not None
                else np.ones(block.size, np.int64))
        far = np.int8(99)

        for gi, g in enumerate(guides):
            best = np.where(ok_f, spacer_mismatches(block, g), far)
            np.minimum(best, np.where(ok_r, spacer_mismatches(rc, g), far), out=best)
            for t in thresholds:
                hit = best <= t
                if not hit.any():
                    continue
                distinct[gi][t] += int(hit.sum())
                sites[gi][t] += int(mult[hit].sum())

        if progress is not None:
            progress(ci + 1, n_chunks)

    return [{"guide_index": i, "distinct": distinct[i], "sites": sites[i]}
            for i in range(len(guides))]


def site_count(universe_codes, universe_mult, query_codes):
    """Both-strand site count in the reference for each query code; 0 when absent."""
    universe_codes = np.asarray(universe_codes, dtype=np.uint64)
    universe_mult = np.asarray(universe_mult, dtype=np.int32)
    q = canonical(query_codes)
    out = np.zeros(q.size, dtype=np.int32)
    if universe_codes.size == 0:
        return out
    idx = np.searchsorted(universe_codes, q)
    idx_clipped = np.clip(idx, 0, universe_codes.size - 1)
    hit = universe_codes[idx_clipped] == q
    out[hit] = universe_mult[idx_clipped[hit]]
    return out
