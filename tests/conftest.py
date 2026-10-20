"""Shared fixtures. Everything here is synthetic, so the suite runs with no genome on disk."""
import gzip

import numpy as np
import pytest


def _wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


@pytest.fixture
def write_fasta(tmp_path):
    """Write records to a FASTA file and return its path. Gzips when the name ends in .gz."""
    def _write(name, records):
        path = tmp_path / name
        body = "".join(f">{rid}\n{_wrap(seq)}\n" for rid, seq in records)
        if name.endswith(".gz"):
            with gzip.open(path, "wt") as fh:
                fh.write(body)
        else:
            path.write_text(body)
        return str(path)
    return _write


@pytest.fixture
def random_sequence():
    """Deterministic random ACGT sequence of a requested length."""
    def _make(length, seed=0):
        rng = np.random.default_rng(seed)
        return "".join(rng.choice(list("ACGT"), size=length))
    return _make
