# Same as compare_panel.py, reading CRISPRitz outputs from the bulge-1 index (timing window run).
"""Panel locus agreement, CRISPRitz v2.9.0 against SWOffinder, hg38 primary chromosomes, 110 guides.

Scope: at most 4 mismatches and at most one bulge of either kind (SWOffinder --maxE 5 --maxM 4
--maxMB 4 --maxB 1; CRISPRitz union of -bDNA 1 -bRNA 0 and -bDNA 0 -bRNA 1 at -mm 4).
Locus key = (guide, chrom, strand, PAM-side forward coordinate), as in compare_chr22.py.
Bulge kind: a gap in the guide row is a DNA bulge, a gap in the genome row is an RNA bulge (both tools).

RNA-only view: CRISPRitz b01 against SWOffinder loci whose reported (best) alignment has no DNA bulge.
SWOffinder keeps one alignment per locus, so its RNA-only set is a lower bound, not a matched set.

Usage: compare_panel.py   (reads panelc/b10_*.targets.txt, panelc/b01_*.targets.txt, ../hg38/panel_keys.tsv)
"""
import collections
import csv
import glob


def kind(target, text):
    return 'DNA' if '-' in target else ('RNA' if '-' in text else 'X')


def crz(setting):
    loci = collections.defaultdict(set)
    rows = 0
    for p in sorted(glob.glob({'b10':'../timing/crz_panel_b10_untimed/*.targets.txt','b01':'../timing/crz_panel/*.targets.txt'}[setting])):
        for r in csv.reader(open(p), delimiter='\t'):
            if not r or r[0].startswith('#'):
                continue
            crrna, dna, ch, pos, st, mm, bs = r[1], r[2], r[3], int(r[4]), r[6], int(r[7]), int(r[8])
            if mm > 4 or bs > 1:
                continue
            rows += 1
            g = crrna.replace('-', '')[:-3] + 'NGG'
            site_len = len(dna.replace('-', ''))
            loci[(g, ch, st, pos + site_len if st == '+' else pos)].add((mm, bs, kind(crrna, dna)))
    return rows, loci


def swo():
    loci = collections.defaultdict(set)
    for r in csv.DictReader(open('../hg38/panel_keys.tsv'), delimiter='\t'):
        mm, bs = int(r['mismatches']), int(r['bulges'])
        if mm > 4 or bs > 1:
            continue
        loci[(r['guide'], r['chrom'], r['strand'], int(r['pam_coord']))].add(
            (mm, bs, kind(r['aligned_target'], r['aligned_text'])))
    return loci


def best(s):
    return min((m + b, b, m) for m, b, _ in s)


def report(name, A, B):
    both, a_only, b_only = set(A) & set(B), set(A) - set(B), set(B) - set(A)
    print(f'[{name}] CRISPRitz loci {len(A)} SWOffinder loci {len(B)} both {len(both)} '
          f'CRISPRitz only {len(a_only)} SWOffinder only {len(b_only)}')
    print(f'[{name}] shared loci with the same best (edits, bulges, mismatches):',
          sum(best(A[k]) == best(B[k]) for k in both), 'of', len(both))
    print(f'[{name}] CRISPRitz-only by best:', collections.Counter(best(A[k]) for k in a_only).most_common(8))
    print(f'[{name}] SWOffinder-only by best:', collections.Counter(best(B[k]) for k in b_only).most_common(8))
    for tag, keys, src in (('crz_only', a_only, A), ('swo_only', b_only, B)):
        with open(f'idx1_panel_{name}_{tag}.tsv', 'w') as fh:
            for k in sorted(keys):
                fh.write('\t'.join(map(str, k)) + '\t' + repr(sorted(src[k])) + '\n')
    per = collections.Counter(k[0] for k in both)
    return per


n10, A10 = crz('b10')
n01, A01 = crz('b01')
A = collections.defaultdict(set)
for d in (A10, A01):
    for k, v in d.items():
        A[k] |= v
B = swo()
print('CRISPRitz rows in scope: b10', n10, 'b01', n01)
report('either', A, B)
B_rna = {k: v for k, v in B.items() if all(x[2] != 'DNA' for x in v)}
report('rna_only', A01, B_rna)
gA = collections.Counter(k[0] for k in A)
gB = collections.Counter(k[0] for k in B)
with open('idx1_panel_per_guide.tsv', 'w') as fh:
    fh.write('guide\tcrispritz_loci\tswoffinder_loci\n')
    for g in sorted(set(gA) | set(gB)):
        fh.write(f'{g}\t{gA[g]}\t{gB[g]}\n')
print('guides with loci: CRISPRitz', len(gA), 'SWOffinder', len(gB),
      'guides where counts differ', sum(gA[g] != gB[g] for g in set(gA) | set(gB)))
