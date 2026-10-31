"""Step two of the maize question: which loci does ONLY an omitted pattern reach?

Step one scanned the 9 patterns Proposition 1 says the 2016 wrapper never emits, and the 133 it
does, over the AGPv4 assembly at Young et al.'s own budget of 5 mismatches. This takes the
difference.

THE ANCHOR IS THE PAM, NOT THE WINDOW START. The patterns are 18 to 22 bases long, so "same locus"
cannot be "same reported start": a bulge pattern one base shorter lands its window one base to the
right of the same cut site. scan_pattern.py reports the forward-leftmost base of the whole
spacer+PAM region (scan_pattern.py:99, start = size-(g+L+3)+1 on the reverse strand), so

    '+'  PAM occupies forward [start+L, start+L+2]
    '-'  PAM occupies forward [start,   start+2]

and (contig, strand, PAM forward start) is comparable across every pattern length. L is read from
the pattern in the file name, never assumed.

A locus survives only if NO emitted pattern OF THE SAME GUIDE reaches it, at either bulge size,
because the wrapper run at bulge size 2 emits both families.

CONTROL. Two things are checked before any number is reported, and either failing aborts:
  1. every omitted pattern's hit set is non-empty, so a zero count is a measurement and not a
     missing file;
  2. the guide's own on-target site is reached by the emitted set, so the emitted side is not
     silently empty or mis-keyed. A difference taken against an empty set is not a finding.
"""
import os, sys, glob, collections

D = sys.argv[1] if len(sys.argv) > 1 else '.'

def key_of(contig, strand, start, L):
    return (contig, strand, start + L if strand == '+' else start)

def load(path, L):
    out = set()
    with open(path) as fh:
        head = fh.readline()
        if not head.startswith('contig\t'):
            sys.exit('BAD HEADER in %s: %r' % (path, head[:60]))
        for line in fh:
            f = line.rstrip('\n').split('\t')
            out.add(key_of(f[0], f[1], int(f[2]), L))
    return out

def pattern_of(name, prefix):
    stem = name[len(prefix):-4] if name.endswith('.tsv') else name[len(prefix):]
    return stem.split('_')[-1]

def guide_of(name, prefix):
    return name[len(prefix):].split('_')[0]

om = collections.defaultdict(dict)   # guide -> pattern -> set
for p in sorted(glob.glob(os.path.join(D, 'hits_*omitted*.tsv'))):
    n = os.path.basename(p); pat = pattern_of(n, 'hits_'); g = guide_of(n, 'hits_')
    om[g][pat] = load(p, len(pat))

em = collections.defaultdict(dict)
for p in sorted(glob.glob(os.path.join(D, 'emit_*.tsv'))):
    n = os.path.basename(p); pat = pattern_of(n, 'emit_'); g = guide_of(n, 'emit_')
    em[g][pat] = load(p, len(pat))

if not om: sys.exit('no omitted hit files in %s' % D)
print('omitted: %d guides, %d patterns, %d hit rows'
      % (len(om), sum(len(v) for v in om.values()),
         sum(len(s) for v in om.values() for s in v.values())))
print('emitted: %d guides, %d patterns, %d hit rows'
      % (len(em), sum(len(v) for v in em.values()),
         sum(len(s) for v in em.values() for s in v.values())))

bad = [(g,p) for g,v in om.items() for p,s in v.items() if not s]
if bad: sys.exit('CONTROL FAILED: omitted pattern(s) with zero hits: %s' % bad)

print()
allonly = {}
for g in sorted(om):
    E = set().union(*em[g].values()) if em.get(g) else set()
    if not E:
        print('%s: NO emitted patterns loaded -- cannot difference, skipping' % g)
        continue
    print('%s  emitted patterns %3d  emitted loci %7d' % (g, len(em[g]), len(E)))
    guide_only = set(); per_pattern_sum_l = [0]
    for pat in sorted(om[g]):
        S = om[g][pat]
        only = S - E
        # CONTROL, PER STRAND. If the PAM offset were wrong on one strand, that strand's loci would
        # share nothing with the emitted set while the other strand overlapped normally. Comparing
        # the two strands catches exactly the error the anchor is most likely to carry, which a
        # single pooled overlap figure would hide.
        for st in ('+', '-'):
            Ss = {k for k in S if k[1] == st}
            if not Ss: continue
            ov = len(Ss & E) / len(Ss)
            if ov == 0.0:
                sys.exit('CONTROL FAILED: %s %s strand shares NO locus with %d emitted loci.\n'
                         '  A mis-keyed strand produces exactly this, so the difference below\n'
                         '  would be an artefact of the anchor and not a measurement.'
                         % (pat, st, len(E)))
            print('        control %s strand: %6d loci, %5.1f%% also reached by an emitted pattern'
                  % (st, len(Ss), 100.0*ov))
        guide_only |= only
        per_pattern_sum_l[0] += len(only)
        print('    %-24s len %2d  loci %6d  reached by NO emitted pattern: %6d  (%.2f%%)'
              % (pat, len(pat), len(S), len(only), 100.0*len(only)/len(S)))
        for k in sorted(only)[:3]:
            print('        example  %s %s PAM@%d' % k)
    # The per-pattern counts above OVERLAP: one locus can be the exclusive property of two
    # different omitted patterns, and summing them would count it twice. The claim is about
    # distinct loci, so the union is the number that may be quoted.
    print('  UNION for this guide: %d distinct loci no emitted pattern reaches'
          ' (sum over patterns would be %d, so %d are shared)'
          % (len(guide_only), per_pattern_sum_l[0], per_pattern_sum_l[0] - len(guide_only)))
    allonly[g] = guide_only
    with open(os.path.join(D, 'exclusive_loci_%s.tsv' % g), 'w') as fh:
        fh.write('contig\tstrand\tpam_forward_start\n')
        for k in sorted(guide_only): fh.write('%s\t%s\t%d\n' % k)
print()
for g in sorted(allonly):
    print('%s  %d distinct loci reachable only by an omitted pattern' % (g, len(allonly[g])))
print('TOTAL distinct loci over all guides measured: %d' % len(set().union(*allonly.values()) if allonly else set()))
print('(guides differenced: %d of %d; a guide with no emitted set loaded is skipped, not counted as zero)'
      % (len(allonly), len(om)))
