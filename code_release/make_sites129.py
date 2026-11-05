#!/usr/bin/env python3
"""Write data/sites129.json from the CHANGE-seq deposit: the capture step of the chain.

This is the producer that was missing from the release. Three scripts read sites129.json
(verify_missed_actives.py, verify_guide_screen.py, validate_terminal_law.py) and nothing
wrote it, so the deposit-to-260 step was auditable by counts and hashes but not re-executable.

The selection is not new. It is the condition already shipped in verify_missed_activity.py
line 26, written there as one line:

    if pos != n-1 or sp[n-1] == sp[n-2]: continue

which is two separate requirements on a single-RNA-bulge measured-active row:

  1. the deleted guide base sits at the terminal spacer position, pos == n-1; and
  2. the guide is at risk under the terminal law, sp[n-1] != sp[n-2].

Requirement 1 alone selects 583 rows over 76 guides. Requirement 2 is what takes those 583
to the 260 over 48 guides in sites129.json; the 323 it removes are the guides the law
rescues, and 323 + 260 = 583 exactly. Dropping requirement 2 is the one way to get a count
that looks like a reproduction failure, so it is spelled out here rather than left in a
conjunction a reader has to split.

Usage, with the deposit on stdin (the same convention as verify_missed_activity.py):

    cat <path-to>/CHANGEseq_CR_Lazzarotto_2020_dataset.csv \
      | python3 code_release/make_sites129.py --out code_release/data/sites129.json

With --check it writes nothing and instead compares what it would have written against the
file already on disk, record for record.
"""
import argparse, csv, hashlib, json, sys

N = 20  # spacer length; the deposit's guides are 20 nt plus a 3 nt PAM


def fingerprint(records):
    """A stable digest of the record SET, independent of file order."""
    keys = sorted('|'.join([r['chrom'], r['sgRNA'], r['aligned_target'],
                            str(r['mm']), '%.6f' % float(r['reads'])]) for r in records)
    if len(set(keys)) != len(keys):
        sys.exit('GATE FAILED: the selected rows are not distinct under the comparison key,\n'
                 '  so a fingerprint match would not establish set identity. Refusing to report one.')
    return hashlib.sha256('\n'.join(keys).encode()).hexdigest()


def select(stream):
    rows = active = bulge = terminal = 0
    out = []
    for r in csv.DictReader(stream):
        rows += 1
        try:
            reads = float((r.get('reads') or '0') or 0)
            nb = int(float((r.get('Align.#Bulges') or '0').strip()))
            nm = int(float((r.get('Align.#Mismatches') or '0').strip()))
        except (TypeError, ValueError):
            continue
        if reads <= 0:
            continue
        active += 1
        if nb < 1:
            continue
        bulge += 1

        g = (r.get('Align.sgRNA') or '').strip().upper()
        t = (r.get('Align.off-target') or '').strip().upper()
        # a single RNA bulge: exactly one gap, and it is in the TARGET, not the guide
        if nb != 1 or t.count('-') != 1 or g.count('-') != 0 or nm > 4:
            continue

        # the deleted guide base's spacer position, counted in real (ungapped) guide bases
        pos = len(g[:t.index('-')].replace('-', ''))
        if pos != N - 1:
            continue                      # requirement 1: terminal position
        terminal += 1

        sp = g.replace('-', '')[:N]
        if sp[N - 1] == sp[N - 2]:
            continue                      # requirement 2: the law's at-risk indicator

        out.append({
            'aligned_target': t,
            'chrom': r.get('chrom'),
            'deleted_position': pos,
            'mm': str(nm),
            'reads': reads,
            'sgRNA': g,
            'spacer': sp,
            'start': r.get('Align.chromStart'),
            'strand': r.get('Align.strand'),
            'target': t.replace('-', ''),
            'window': r.get('SiteWindow'),
        })
    if rows == 0:
        sys.exit("GATE FAILED: no rows read from stdin, so nothing was selected. This is an\n"
                 "  invocation error, not a finding. Pipe the CHANGE-seq deposit in, e.g.\n"
                 "  cat <path-to>/"
                 "CHANGEseq_CR_Lazzarotto_2020_dataset.csv \\\n"
                 "    | python3 code_release/make_sites129.py --out code_release/data/sites129.json")
    return out, dict(rows=rows, active=active, bulge=bulge, terminal=terminal)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', help='path to write; omit with --check')
    ap.add_argument('--check', metavar='PATH',
                    help='compare against an existing file instead of writing')
    a = ap.parse_args()
    if not a.out and not a.check:
        ap.error('give --out or --check')

    recs, stage = select(sys.stdin)
    guides = len({r['spacer'] for r in recs})
    print('deposit rows %(rows)d -> measured active %(active)d -> bulge-bearing %(bulge)d' % stage)
    print('  -> single RNA bulge at the terminal position within 4 mm: %d  (requirement 1)'
          % stage['terminal'])
    print('  -> of those, at risk under the terminal law: %d over %d guides  (requirement 2)'
          % (len(recs), guides))
    print('  -> rescued by the law and therefore excluded: %d' % (stage['terminal'] - len(recs)))
    fp = fingerprint(recs)
    print('  record-set SHA-256: %s' % fp)

    if a.check:
        held = json.load(open(a.check))
        hfp = fingerprint(held)
        same = hfp == fp and len(held) == len(recs)
        print('\nheld file %s: %d records, SHA-256 %s' % (a.check, len(held), hfp))
        print('SET IDENTICAL: %s' % same)
        sys.exit(0 if same else 'MISMATCH: the selection does not reproduce the held file.')

    with open(a.out, 'w') as fh:
        json.dump(recs, fh, indent=1, sort_keys=True)
    print('\nwrote %d records to %s' % (len(recs), a.out))


if __name__ == '__main__':
    main()
