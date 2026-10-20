"""Two-bit coding of SpCas9 target windows, and the one distinction that governs correctness.

A SpCas9 target is a 23-mer: a 20-base spacer followed by an NGG PAM. We pack it into the low
46 bits of a uint64, two bits per base, most significant base first. So base 1 of the spacer
occupies bits 45-44 and base 23 (the second G of the PAM) occupies bits 1-0.

**Oriented versus canonical, which is the whole game.**

A guide binds double-stranded DNA, so one 23-mer of sequence corresponds to a binding site on
whichever strand presents it as ``[spacer]NGG``. Two different operations need two different
representations of that site, and conflating them silently corrupts a comparison:

- **Membership** in a reference genome is strand-free. A site is present if the reference
  contains it on either strand, so membership is tested on the *canonical* code,
  ``min(code, revcomp(code))``.
- **Mismatch distance** to a guide spacer is strand-specific. It must be computed on the
  *oriented* code, the one that actually reads ``[spacer]NGG``.

About 58% of canonical codes are the reverse-complement of their oriented form. A canonical
code therefore begins with ``CC`` roughly half the time, which is a compositional constraint
that no real spacer has. Using canonical codes for distance imposes that constraint on whichever
arm of a comparison you canonicalised, and not on the other. This package keeps the two strictly
separate: :func:`oriented_sites` returns oriented codes, and canonicalisation happens only inside
the membership test.
"""
from __future__ import annotations

import gzip

import numpy as np

#: Length of an SpCas9 target window: 20-base spacer plus NGG PAM.
K = 23
#: Bases of PAM at the 3' end of the window.
PAM_LEN = 3
#: Bases of spacer at the 5' end of the window.
SPACER_LEN = K - PAM_LEN
BITS = 2 * K

_M46 = np.uint64((1 << BITS) - 1)
_SB = np.uint64(64 - BITS)
_R2 = np.uint64(0x3333333333333333)
_R4 = np.uint64(0x0F0F0F0F0F0F0F0F)
_R8 = np.uint64(0x00FF00FF00FF00FF)
_R16 = np.uint64(0x0000FFFF0000FFFF)
_S2, _S4, _S8, _S16, _S32 = (np.uint64(x) for x in (2, 4, 8, 16, 32))

#: Base-to-value table. A=0, C=1, G=2, T=3; everything else, N included, is 255 and invalid.
CODE_TABLE = np.full(256, 255, np.uint8)
for _b, _v in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    CODE_TABLE[_b] = _v

_A, _C, _G, _T = 0, 1, 2, 3


def revcomp64(x):
    """Reverse-complement packed 23-mer codes, elementwise.

    Complement is a bitwise negation within 46 bits because A/T and C/G are 0/3 and 1/2.
    Reversal is a sequence of swaps over 2-bit groups. Accepts a scalar or an array of uint64.
    """
    x = np.asarray(x, dtype=np.uint64)
    x = x ^ _M46
    x = ((x & _R2) << _S2) | ((x >> _S2) & _R2)
    x = ((x & _R4) << _S4) | ((x >> _S4) & _R4)
    x = ((x & _R8) << _S8) | ((x >> _S8) & _R8)
    x = ((x & _R16) << _S16) | ((x >> _S16) & _R16)
    x = (x << _S32) | (x >> _S32)
    return x >> _SB


def canonical(codes):
    """The strand-free form of a code, for membership testing only.

    Never use the result for mismatch distance; see this module's docstring.
    """
    codes = np.asarray(codes, dtype=np.uint64)
    return np.minimum(codes, revcomp64(codes))


def encode_window(seq: str):
    """Pack one 23-mer string into a code, or return ``None`` if it is not encodable.

    Returns ``None`` for the wrong length or any non-ACGT character, rather than raising, so
    callers can skip malformed rows in a panel file without special-casing.
    """
    vals = CODE_TABLE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    if vals.size != K or (vals > 3).any():
        return None
    code = np.uint64(0)
    for v in vals:
        code = (code << np.uint64(2)) | np.uint64(int(v))
    return int(code)


def encode_spacer(spacer: str):
    """Pack a bare 20-base spacer into a window code, or ``None`` if unencodable.

    The spacer is placed in the spacer bits and the PAM bits are left zero. That is exact rather
    than a placeholder: :func:`scotoma.distance.spacer_mismatches` discards the PAM bits before
    comparing, so a guide's own PAM never affects a distance. This exists so a panel can be
    given as 20-mers, which is how guides are usually published, without inventing an NGG.
    """
    vals = CODE_TABLE[np.frombuffer(spacer.encode(), dtype=np.uint8)]
    if vals.size != SPACER_LEN or (vals > 3).any():
        return None
    code = np.uint64(0)
    for v in vals:
        code = (code << np.uint64(2)) | np.uint64(int(v))
    return int(code << np.uint64(2 * PAM_LEN))


def decode_window(code: int) -> str:
    """Unpack a code back to its 23-mer string. Inverse of :func:`encode_window`."""
    code = int(code)
    out = []
    for i in range(K - 1, -1, -1):
        out.append("ACGT"[(code >> (2 * i)) & 0b11])
    return "".join(out)


def spacer_of(code: int) -> str:
    """The 20-base spacer of a packed window, as a string."""
    return decode_window(code)[:SPACER_LEN]


def spacer_gc(codes):
    """Number of G or C bases in the 20-base spacer of each code, as int8 in ``0 .. 20``.

    Computed on the packed representation without decoding. With A=00, C=01, G=10 and T=11, a
    base is G or C exactly when its two bits differ, so XORing each 2-bit group's high bit onto
    its low bit and population-counting the low bits gives the GC count directly. The PAM is
    excluded, matching :func:`scotoma.distance.spacer_mismatches`, so the statistic describes the
    part of the window a guide recognises.
    """
    codes = np.asarray(codes, dtype=np.uint64)
    x = codes >> np.uint64(2 * PAM_LEN)
    y = ((x >> np.uint64(1)) ^ x) & np.uint64(0x5555555555555555)
    c = y - ((y >> np.uint64(1)) & np.uint64(0x5555555555555555))
    c = (c & np.uint64(0x3333333333333333)) + ((c >> np.uint64(2)) & np.uint64(0x3333333333333333))
    c = (c + (c >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return ((c * np.uint64(0x0101010101010101)) >> np.uint64(56)).astype(np.int8)


def shuffle_spacer(code: int, rng) -> int:
    """Permute a window's spacer bases, preserving its exact base composition.

    The PAM bits are left untouched; only the 20 spacer positions are reordered. Composition is
    preserved exactly, so GC content and mononucleotide counts are identical to the original while
    base order is destroyed.

    **What this controls.** A composition-preserving shuffled panel is the negative control for a
    matched-arm density comparison. If the candidate-density deficit between two arms is a property
    of the *arms*, shuffling the guides must leave it intact. If shuffling abolishes it, the deficit
    is specific to the real guide sequences, which matters here because a published guide panel was
    designed and measured against the reference: real guides could be closer to reference-present
    sequence by provenance rather than by biology.
    """
    bases = list(decode_window(int(code))[:SPACER_LEN])
    order = rng.permutation(len(bases))
    shuffled = "".join(bases[i] for i in order)
    pam_bits = int(code) & ((1 << (2 * PAM_LEN)) - 1)
    packed = encode_spacer(shuffled)
    if packed is None:  # pragma: no cover - a permutation of ACGT is always encodable
        raise ValueError(f"could not encode shuffled spacer {shuffled!r}")
    return int(packed) | pam_bits


def shuffle_panel(guide_codes, seed: int = 0):
    """Composition-preserving shuffle of every guide in a panel. Returns a new uint64 array."""
    rng = np.random.default_rng(seed + 7000)
    return np.array([shuffle_spacer(int(g), rng) for g in np.asarray(guide_codes).tolist()],
                    dtype=np.uint64)


def window_complexity(codes, word: int = 3, chunk: int = 2_000_000):
    """Distinct ``word``-mer count in each 23-mer window: a linguistic-complexity statistic.

    Returns an int8 in ``1 .. K - word + 1``. For the default 3-mers that is 1 to 21: a perfect
    homopolymer gives 1, and a complex window approaches 21.

    **Why this exists.** Matching two arms on GC content does not control for higher-order
    repetitiveness. Sequence absent from a reference is heavily enriched for alpha-satellite,
    HSat2 and HSat3 arrays, rDNA and simple tandem repeats, all of which are periodic and
    low-entropy. A complex 20-base spacer finds fewer near-matches in low-entropy sequence for
    reasons that have nothing to do with whether a reference can represent that sequence. So a
    deficit in candidate density between reference-absent and reference-present arms can be
    produced by complexity alone, and a GC-matched comparison cannot tell the two apart. This
    statistic exists to be matched on.

    Computed by extracting all ``K - word + 1`` overlapping words per window, sorting each
    window's words, and counting distinct runs. Chunked so peak memory stays near
    ``chunk * (K - word + 1)`` bytes.
    """
    codes = np.asarray(codes, dtype=np.uint64)
    n_words = K - word + 1
    if n_words < 1:
        raise ValueError(f"word must be at most K={K}")
    width = np.uint64(2 * word)
    mask = np.uint64((1 << int(width)) - 1)
    out = np.empty(codes.size, np.int8)

    for start in range(0, codes.size, chunk):
        stop = min(start + chunk, codes.size)
        block = codes[start:stop]
        words = np.empty((block.size, n_words), np.uint16)
        # Word j occupies bits [2*(K-word-j) .. ) counting from the low end, because base 1 is
        # the most significant base.
        for j in range(n_words):
            shift = np.uint64(2 * (K - word - j))
            words[:, j] = ((block >> shift) & mask).astype(np.uint16)
        words.sort(axis=1)
        distinct = np.ones(block.size, np.int8)
        distinct += (words[:, 1:] != words[:, :-1]).sum(axis=1).astype(np.int8)
        out[start:stop] = distinct
        del words
    return out


def iter_fasta(path):
    """Yield ``(record_name, bases)`` for a FASTA file, gzipped or not.

    ``bases`` is a uint8 numpy array of the raw ASCII, suitable for indexing
    :data:`CODE_TABLE`. Records are concatenated across lines; whitespace is stripped.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None and chunks:
                    yield name, np.frombuffer(b"".join(chunks), np.uint8)
                name = line[1:].split()[0].decode()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None and chunks:
        yield name, np.frombuffer(b"".join(chunks), np.uint8)


#: The PAM scopes this package can extract. ``"NGG"`` is the canonical SpCas9 PAM.
#: ``"NGG+NAG"`` additionally admits the weaker NAG PAM, which off-target screens often include
#: and which materially changes the number of sites: NGG alone is 1 window in 16 per strand,
#: NGG or NAG is 1 in 8, so a site count under the wider scope is roughly twice as large.
#:
#: **A count under one scope is not comparable to a count under the other.** Mixing them is a
#: denominator error, and this package records the scope in every universe it builds so that a
#: later comparison can be refused rather than silently misread.
PAM_SCOPES = ("NGG", "NGG+NAG")
DEFAULT_PAM = "NGG"


def oriented_sites(bases, pam: str = DEFAULT_PAM):
    """Extract every SpCas9 binding site in a sequence as an **oriented** 23-mer code.

    A site exists wherever a 23-base window reads ``[spacer][PAM]`` on either strand. Under the
    default ``NGG`` scope that is, on the forward strand, a window whose last two bases are
    ``GG``; on the reverse strand it appears in the forward sequence as a window whose *first*
    two bases are ``CC``, and its oriented code is the reverse complement of that window. Under
    ``NGG+NAG`` the 22nd base may also be ``A`` (and correspondingly ``T`` at the second position
    of a reverse-strand window).

    Windows containing any non-ACGT base are dropped. A window that satisfies the PAM rule on
    both strands yields two codes, because it genuinely is two binding sites.

    Parameters
    ----------
    bases : numpy.ndarray
        uint8 ASCII, as yielded by :func:`iter_fasta`.
    pam : str
        One of :data:`PAM_SCOPES`.

    Returns
    -------
    numpy.ndarray
        uint64 oriented codes, unsorted, with multiplicity. Empty if the sequence is shorter
        than :data:`K` or contains no valid site.
    """
    if pam not in PAM_SCOPES:
        raise ValueError(f"pam must be one of {PAM_SCOPES}, got {pam!r}")
    bases = np.asarray(bases, dtype=np.uint8)
    vals = CODE_TABLE[bases]
    bad = vals == 255
    v = np.where(bad, 0, vals).astype(np.uint64)
    n = v.size - K + 1
    if n <= 0:
        return np.empty(0, np.uint64)

    codes = np.zeros(n, np.uint64)
    for j in range(K):
        codes = (codes << np.uint64(2)) | v[j:j + n]

    # A window is usable only if it contains no invalid base. Prefix sums over `bad` give that
    # test for every window at once.
    bc = np.cumsum(np.concatenate(([0], bad.astype(np.int64))))
    ok = (bc[K:K + n] - bc[0:n]) == 0

    b22, b23 = v[K - 2:K - 2 + n], v[K - 1:K - 1 + n]
    r1, r2 = v[0:n], v[1:1 + n]
    if pam == "NGG":
        fwd_pam = ok & (b22 == _G) & (b23 == _G)
        rev_pam = ok & (r1 == _C) & (r2 == _C)
    else:
        # NGG or NAG. On the reverse strand the complement of [G or A] at position 22 is
        # [C or T] at position 2, and the complement of G at 23 is C at 1.
        fwd_pam = ok & ((b22 == _G) | (b22 == _A)) & (b23 == _G)
        rev_pam = ok & (r1 == _C) & ((r2 == _C) | (r2 == _T))

    parts = [codes[np.flatnonzero(fwd_pam)]]
    rev = codes[np.flatnonzero(rev_pam)]
    if rev.size:
        parts.append(revcomp64(rev))
    return np.concatenate(parts) if parts else np.empty(0, np.uint64)
