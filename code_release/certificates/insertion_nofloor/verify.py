"""Read-only exact replay of insertion covers and finite lower certificates."""

from collections import Counter
from itertools import combinations, product
from math import comb
from pathlib import Path
import gzip
import csv
import json
import sys

ROOT = Path(__file__).resolve().parent
ALPHABET = "ACGT"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def word_list(n):
    return ["".join(x) for x in product(ALPHABET, repeat=n)]


def distance(a, b):
    require(len(a) == len(b), "unequal word lengths")
    return sum(a[i] != b[i] for i in range(len(a)))


def deletion_oracle(x, guide):
    """Distance to the insertion family, computed without building insertions."""
    return min(distance(x[:j] + x[j + 1:], guide) for j in range(len(x)))


def coverage_problem(guide, m):
    require(m >= 0 and all(c in ALPHABET for c in guide), "invalid problem")
    space = word_list(len(guide) + 1)
    lib = [x for x in space if deletion_oracle(x, guide) == 0]
    points = [x for x in space if deletion_oracle(x, guide) <= m]
    masks = {c: sum(1 << j for j, x in enumerate(points) if distance(c, x) <= m)
             for c in space}
    return space, lib, points, masks


def set_bits(mask):
    while mask:
        b = mask & -mask
        yield b.bit_length() - 1
        mask -= b


def check_cover(data):
    space, lib, points, masks = coverage_problem(data["guide"], data["m"])
    cover = data["cover"]
    require(len(cover) == len(set(cover)), "duplicate cover centres")
    require(all(c in masks for c in cover), "invalid cover word")
    covered = 0
    for c in cover:
        covered |= masks[c]
    if covered != (1 << len(points)) - 1:
        missing = next(x for i, x in enumerate(points) if not covered >> i & 1)
        raise ValueError("uncovered target " + missing)
    return space, lib, points, masks


def check_certificate(data):
    space, lib, points, masks = check_cover(data)
    require(data["family"] == lib, "insertion family disagrees with deletion oracle")
    require(len(lib) == 3 * len(data["guide"]) + 4, "wrong count identity")
    require(data["space_size"] == len(space), "wrong ambient space size")
    require(data["universe_size"] == len(points), "wrong target size")
    require(data["K"] == len(data["cover"]), "wrong claimed cover size")
    require(data["infeasible_k"] == data["K"] - 1, "lower bound does not match upper")
    reps = data["retained_centres"]
    require(len(reps) == len(set(reps)), "duplicate retained centre")
    require(all(c in masks for c in reps), "invalid retained centre")
    reduced = [masks[c] for c in reps]
    dom = data["dominators"]
    require(len(dom) == len(space), "not all four-letter centres accounted for")
    for c, j in zip(space, dom):
        if j == -1:
            require(not masks[c], "nonempty centre discarded without a dominator")
        else:
            require(isinstance(j, int) and 0 <= j < len(reps), "invalid dominator index")
            require(masks[c] & ~reduced[j] == 0, "discarded centre has no covering superset")
    point_options = [sum(1 << j for j, mask in enumerate(reduced) if mask >> p & 1)
                     for p in range(len(points))]
    proof = data["lower_proof"]
    nodes = proof["nodes"]
    visited, active = {}, set()
    leaf_counts = Counter()

    def visit(node_id, rem, avail, slots):
        require(isinstance(node_id, int) and 0 <= node_id < len(nodes), "bad proof node")
        state = rem, avail, slots
        if node_id in visited:
            require(visited[node_id] == state, "proof node reused for a different state")
            return
        require(node_id not in active, "cyclic proof")
        active.add(node_id)
        require(rem != 0, "proof incorrectly rejects a complete cover")
        node = nodes[node_id]
        kind = node[0]
        leaf_counts[kind] += 1
        if kind == "zero":
            require(slots <= 0, "false zero-budget leaf")
        elif kind == "volume":
            require(slots >= 0, "invalid volume budget")
            gains = sorted(((reduced[i] & rem).bit_count() for i in set_bits(avail)), reverse=True)
            require(sum(gains[:slots]) < rem.bit_count(), "false volume bound")
        elif kind == "empty":
            p = node[1]
            require(0 <= p < len(points) and rem >> p & 1, "empty leaf uses a covered target")
            require(point_options[p] & avail == 0, "target still has an available centre")
        elif kind == "packing":
            witness = node[1]
            require(len(witness) == len(set(witness)) and len(witness) > slots,
                    "packing is too small or has duplicate points")
            occupied = 0
            for p in witness:
                require(0 <= p < len(points) and rem >> p & 1, "packing point is not uncovered")
                options = point_options[p] & avail
                require(options != 0, "packing point has no candidate")
                require(occupied & options == 0, "one centre covers two packing points")
                occupied |= options
        elif kind == "branch":
            p, children = node[1:]
            require(0 <= p < len(points) and rem >> p & 1, "branch target is already covered")
            choices = point_options[p] & avail
            branch_indices = [item[0] for item in children]
            require(len(branch_indices) == len(set(branch_indices)), "duplicate proof branch")
            require(set(branch_indices) == set(set_bits(choices)), "proof omits an available branch")
            require(slots > 0 and choices, "invalid branch budget or pivot")
            for i, child in children:
                avail &= ~(1 << i)
                visit(child, rem & ~reduced[i], avail, slots - 1)
        else:
            raise ValueError("unknown proof node type")
        active.remove(node_id)
        visited[node_id] = state

    visit(proof["root"], (1 << len(points)) - 1,
          (1 << len(reps)) - 1, data["infeasible_k"])
    require(len(visited) == len(nodes), "certificate contains unreachable proof nodes")
    return {"status": "VERIFIED-COMPUTATIONALLY", "guide": data["guide"], "m": data["m"],
            "K": data["K"], "targets": len(points), "all_centres": len(space),
            "proof_nodes": len(visited), "node_types": dict(leaf_counts)}


def read_certificate(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def count_checks(max_n=7):
    checked = 0
    coordinate_checks = 0
    for n in range(max_n + 1):
        for g in word_list(n):
            lib = {g[:i] + a + g[i:] for i in range(n + 1) for a in ALPHABET}
            require(len(lib) == 3 * n + 4, "count identity failed")
            patterns = {g[:i] + "*" + g[i:] for i in range(n + 1)}
            require(len(patterns) == n + 1, "wildcard count failed")
            if n:
                nondegenerate = {g[:i] + a + g[i:] for i in range(1, n + 1) for a in ALPHABET}
                require(len(nondegenerate) == 3 * n + 1, "nondegenerate count failed")
            for j in range(n + 1):
                require({x[j] for x in lib} == set(ALPHABET), "coordinate saturation failed")
                coordinate_checks += 1
            checked += 1
    return {"status": "VERIFIED-COMPUTATIONALLY", "n_range": [0, max_n],
            "guides": checked, "saturated_coordinates": coordinate_checks}


def runs(g):
    lengths = []
    for i, a in enumerate(g):
        if not i or a != g[i - 1]:
            lengths.append(1)
        else:
            lengths[-1] += 1
    return tuple(lengths)


def isometry_checks(rows):
    refs = [r["guide"] for r in rows if r["m"] == 1]
    mapping_counts = Counter()
    checked = 0
    for n in range(1, 5):
        for guide in word_list(n):
            candidates = [g for g in refs if runs(g) in (runs(guide), runs(guide)[::-1])]
            require(candidates, "table omits a guide run class")
            ref = candidates[0]
            source = guide if runs(guide) == runs(ref) else guide[::-1]
            maps = []
            for j in range(n + 1):
                sb = [source[t] for t in (j - 1, j) if 0 <= t < n]
                rb = [ref[t] for t in (j - 1, j) if 0 <= t < n]
                perm = dict(zip(sb, rb))
                unused_s = [a for a in ALPHABET if a not in perm]
                unused_r = [a for a in ALPHABET if a not in perm.values()]
                perm.update(zip(unused_s, unused_r))
                require(len(perm) == len(set(perm.values())) == 4, "not a coordinate permutation")
                maps.append(perm)
            lib_source = {source[:i] + a + source[i:] for i in range(n + 1) for a in ALPHABET}
            lib_ref = {ref[:i] + a + ref[i:] for i in range(n + 1) for a in ALPHABET}
            image = {"".join(maps[j][a] for j, a in enumerate(x)) for x in lib_source}
            require(image == lib_ref, "insertion-family isometry failed")
            checked += 1
            mapping_counts[ref] += 1
    require(checked == 340, "not every guide of lengths one through four was checked")
    return {"status": "VERIFIED-COMPUTATIONALLY", "guides": checked,
            "representative_multiplicities": dict(mapping_counts)}


def long_witness_check():
    data = json.loads((ROOT / "long_witness.json").read_text())
    n, m = len(data["guide"]), data["m"]
    require(n == 20 and data["guide"] == "A" * n and m == 4,
            "wrong length-20 homopolymer problem")
    N = n + 1
    base = "A" * N
    chosen = data["chosen_positions_zero_based"]
    require(len(chosen) == len(set(chosen)) == N - m, "wrong chosen-position count")
    require(all(1 <= p < N for p in chosen), "chosen positions include degenerate slot")
    concrete = {base}
    for p in chosen:
        for a in "CGT":
            concrete.add(base[:p] + a + base[p + 1:])
    patterns = {base[:p] + "*" + base[p + 1:] for p in chosen}
    require(set(data["concrete_cover"]) == concrete, "wrong long concrete cover words")
    require(set(data["wildcard_cover"]) == patterns, "wrong long wildcard cover words")
    require(len(data["concrete_cover"]) == len(concrete), "duplicate long concrete centre")
    require(len(data["wildcard_cover"]) == len(patterns), "duplicate long wildcard pattern")
    family = {data["guide"][:i] + a + data["guide"][i:]
              for i in range(N) for a in ALPHABET}
    positive_family = {data["guide"][:i] + a + data["guide"][i:]
                       for i in range(1, N) for a in ALPHABET}
    require(concrete <= positive_family <= family, "long cover is not sound by containment")
    for name, expected in (("long_concrete_cover.txt", concrete),
                           ("long_concrete_family.txt", family),
                           ("long_wildcard_cover.txt", patterns)):
        lines = (ROOT / name).read_text().splitlines()
        require(len(lines) == len(expected) and set(lines) == expected,
                "written long witness differs from certificate: " + name)
    for name, expected in (("concrete_family_size", len(family)),
                           ("nondegenerate_concrete_family_size", len(positive_family)),
                           ("wildcard_family_size", N),
                           ("nondegenerate_wildcard_family_size", n)):
        require(data[name] == expected, "wrong long witness metadata: " + name)
    selected_mask = sum(1 << p for p in chosen)
    full_targets = 0
    production_targets = 0
    concrete_accepted = 0
    wildcard_accepted = 0
    powers = [3 ** w for w in range(N + 1)]
    for support in range(1 << N):
        weight = support.bit_count()
        first = support & 1
        meets = bool(support & selected_mask)
        # Independently derived support-class distance formulas. Their exact
        # translation to the written words is checked above and in small spaces.
        d_family = max(0, weight - 1)
        d_production = first + max(0, weight - first - 1)
        d_concrete = weight - int(meets)
        d_pattern = weight - int(meets)
        target = d_family <= m
        production_target = d_production <= m
        c_accept, p_accept = d_concrete <= m, d_pattern <= m
        require(target == production_target == c_accept == p_accept,
                "support class is missed or spuriously accepted")
        full_targets += powers[weight] * target
        production_targets += powers[weight] * production_target
        concrete_accepted += powers[weight] * c_accept
        wildcard_accepted += powers[weight] * p_accept
    require(full_targets == sum(comb(N, t) * 3 ** t for t in range(m + 2)),
            "incorrect Hamming-ball volume")
    require(data["target_size"] == full_targets == 5467456,
            "wrong long target count")
    # Compare the support formulas with literal word membership, without a quotient,
    # at every nontrivial small homopolymer instance through guide length four.
    literal_targets = 0
    for small_n in range(1, 5):
        small_N = small_n + 1
        for small_m in range(1, small_n + 1):
            small_base = "A" * small_N
            pos = list(range(1, small_N - small_m + 1))
            cc = [small_base] + [small_base[:p] + a + small_base[p + 1:]
                                for p in pos for a in "CGT"]
            pp = [small_base[:p] + "*" + small_base[p + 1:] for p in pos]
            for x in word_list(small_N):
                weight = sum(a != "A" for a in x)
                actual_c = min(distance(x, c) for c in cc) <= small_m
                actual_p = any(sum(a != b for a, b in zip(x, p) if b != "*") <= small_m for p in pp)
                require(actual_c == actual_p == (weight <= small_m + 1), "literal long-rule calibration failed")
                require((deletion_oracle(x, "A" * small_n) <= small_m) == actual_c,
                        "homopolymer alignment oracle failed")
                literal_targets += 1
    return {"status": "VERIFIED-COMPUTATIONALLY", "n": n, "m": m,
            "concrete_family_size": 3 * n + 4, "nondegenerate_concrete_size": 3 * n + 1,
            "concrete_cover_size": len(concrete), "wildcard_family_size": N,
            "nondegenerate_wildcard_size": n, "wildcard_cover_size": len(patterns),
            "ambient_words": 4 ** N, "support_classes_checked": 2 ** N,
            "target_words": full_targets, "small_literal_checks": literal_targets}


def smallest_brute_force():
    """Literal coverage and all ambient triples, separate from the branching DAG."""
    space = word_list(2)
    family = (ROOT / "smallest_family.txt").read_text().splitlines()
    cover = (ROOT / "smallest_cover.txt").read_text().splitlines()
    require(set(family) == {"AA", "AC", "AG", "AT", "CA", "GA", "TA"}
            and len(family) == 7, "wrong smallest family")
    require(cover == ["AA", "AC", "AG", "AT"], "wrong smallest cover")
    masks = {c: sum(1 << i for i, x in enumerate(space) if distance(c, x) <= 1)
             for c in space}
    full = (1 << len(space)) - 1
    require(all(min(distance(x, s) for s in family) <= 1 for x in space),
            "smallest target union is not the full space")
    require(all(min(distance(x, c) for c in cover) <= 1 for x in space),
            "smallest cover misses a target")
    triples = 0
    for a, b, c in combinations(space, 3):
        require(masks[a] | masks[b] | masks[c] != full, "three centres suffice")
        triples += 1
    require(triples == comb(16, 3) == 560, "incomplete triple exhaustion")
    with (ROOT / "smallest_membership.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require([row["target"] for row in rows] == space, "membership table is incomplete")
    for row in rows:
        x = row["target"]
        df = min(distance(x, s) for s in family)
        dc = min(distance(x, c) for c in cover)
        require([int(row[k]) for k in ("distance_to_family", "distance_to_cover",
                                      "in_union", "covered")] == [df, dc, 1, 1],
                "incorrect literal membership table")
    require(all(sum(a != b for a, b in zip(x, "A*") if b != "*") <= 1
                for x in space), "smallest wildcard cover misses a target")
    return {"status": "VERIFIED-COMPUTATIONALLY", "targets": 16,
            "ambient_triples_rejected": triples, "K_word": 4, "K_1star": 1}


def run_extension_check(data=None):
    """Exhaust the removed ball in the additional four-letter guide witness."""
    if data is None:
        data = json.loads((ROOT / "run_extension_witness.json").read_text())
    guide, m = data["guide"], data["m"]
    removed, words = data["removed_centre"], data["cover"]
    cover = set(words)
    require(guide == "A" * 17 + "CGT" and m == 4 and removed == "A" * 18 + "CGT",
            "wrong run-extension problem")
    family = {guide[:i] + a + guide[i:] for i in range(21) for a in ALPHABET}
    require(len(family) == 64 and cover == family - {removed} and len(words) == 63,
            "run-extension cover is not the specified family minus one word")
    written = (ROOT / "run_extension_cover.txt").read_text().splitlines()
    require(len(written) == 63 and set(written) == cover, "written run-extension cover differs")
    neighbours = {(p, a): removed[:p] + a + removed[p + 1:]
                  for p in range(18) for a in "CGT"}
    require(set(neighbours.values()) <= cover, "missing run-extension neighbour")
    count = 0
    for changes in range(m + 1):
        for positions in combinations(range(21), changes):
            alternatives = [[a for a in ALPHABET if a != removed[p]] for p in positions]
            for letters in product(*alternatives):
                x = list(removed)
                for p, a in zip(positions, letters):
                    x[p] = a
                target = "".join(x)
                require(distance(target, removed) == changes, "wrong removed-ball member")
                p = next((p for p in positions if p < 18), None)
                chosen = neighbours[(p, target[p])] if p is not None else neighbours[(0, "C")]
                require(chosen in cover and distance(target, chosen) <= m,
                        "removed-ball target has no retained cover")
                count += 1
    require(count == sum(comb(21, j) * 3 ** j for j in range(5)) == 522649,
            "removed ball was not exhausted")
    return {"status": "VERIFIED-COMPUTATIONALLY", "family_size": 64,
            "cover_size": 63, "removed_ball_targets_exhausted": count}


def expect_failure(name, function):
    try:
        function()
    except (ValueError, KeyError, IndexError, TypeError) as error:
        return {"control": name, "status": "VERIFIED-COMPUTATIONALLY",
                "rejected": True, "reason": str(error)}
    raise ValueError("negative control incorrectly accepted: " + name)


def controls():
    good = read_certificate(ROOT / "certificates/A_m1.json.gz")
    results = []
    bad_cover = {"guide": "A", "m": 1, "cover": ["AC", "AG", "AT"]}
    results.append(expect_failure("remove AA from four-centre cover", lambda: check_cover(bad_cover)))
    bad_radius = {"guide": "A", "m": 0, "cover": ["AA", "AC", "AG", "AT"]}
    results.append(expect_failure("reuse positive-budget cover at zero", lambda: check_cover(bad_radius)))
    incomplete = json.loads(json.dumps(good))
    incomplete["dominators"] = incomplete["dominators"][:-1]
    results.append(expect_failure("omit an ambient centre", lambda: check_certificate(incomplete)))
    bad_branch = json.loads(json.dumps(good))
    root = bad_branch["lower_proof"]["root"]
    require(bad_branch["lower_proof"]["nodes"][root][0] == "branch", "unexpected control root")
    bad_branch["lower_proof"]["nodes"][root][2].pop()
    results.append(expect_failure("omit a lower-proof branch", lambda: check_certificate(bad_branch)))
    bad_leaf = json.loads(json.dumps(good))
    bad_leaf["lower_proof"]["nodes"][root] = ["volume"]
    results.append(expect_failure("plant false volume inequality", lambda: check_certificate(bad_leaf)))
    results.append(expect_failure("replace 3n+4 by 3n+3",
                                  lambda: require(len(coverage_problem("A", 0)[1]) == 6,
                                                  "seven insertion words, not six")))
    # An invalid seven-unit packing at radius one must be rejected as well.
    lib = coverage_problem("A", 0)[1]
    results.append(expect_failure("assign unit weight to the seven library words",
                                  lambda: require(max(sum(distance(c, x) <= 1 for x in lib)
                                                      for c in word_list(2)) <= 1,
                                                  "a ball collects seven units")))
    return results


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "certificate":
        print(json.dumps(check_certificate(read_certificate(Path(sys.argv[2]))), sort_keys=True))
        return
    if len(sys.argv) == 3 and sys.argv[1] == "cover":
        data = json.loads(Path(sys.argv[2]).read_text())
        check_cover(data)
        print("VERIFIED-COMPUTATIONALLY: cover accepted")
        return
    require(sys.argv[1:] in ([], ["verify"]),
            "usage: verify.py [verify | certificate FILE | cover FILE]")
    rows = json.loads((ROOT / "table.json").read_text())
    require(len(rows) == 14 and len({(r["guide"], r["m"]) for r in rows}) == 14,
            "expected 14 distinct exact cases")
    checked = []
    for row in rows:
        data = read_certificate(ROOT / row["certificate"])
        result = check_certificate(data)
        for key in ("guide", "m", "K"):
            require(result[key] == row[key], "table disagrees with certificate: " + key)
        for key, value in (("n", len(data["guide"])), ("family_size", len(data["family"])),
                           ("wildcard_family_size", len(data["guide"]) + 1),
                           ("universe_size", result["targets"]),
                           ("space_size", result["all_centres"]),
                           ("proof_nodes", result["proof_nodes"]),
                           ("retained_centres", len(data["retained_centres"])),
                           ("cover", data["cover"])):
            require(row[key] == value, "table metadata disagrees: " + key)
        checked.append(result)
        print("VERIFIED-COMPUTATIONALLY: guide={guide} m={m} K={K} targets={targets} "
              "ambient_centres={all_centres} proof_nodes={proof_nodes}".format(**result), flush=True)
    with (ROOT / "exact_table.csv").open(newline="") as handle:
        csv_rows = list(csv.DictReader(handle))
    require(len(csv_rows) == 14, "CSV table has wrong number of rows")
    by_case = {(r["guide"], r["m"]): r for r in rows}
    for row in csv_rows:
        ref = by_case[(row["guide"], int(row["m"]))]
        for key in ("n", "m", "K", "family_size", "wildcard_family_size",
                    "universe_size", "space_size", "proof_nodes"):
            require(int(row[key]) == ref[key], "CSV table disagrees: " + key)
        require(row["certificate"] == ref["certificate"], "CSV certificate differs")
    checks = [("count_checks", count_checks), ("isometries", lambda: isometry_checks(rows)),
              ("smallest_brute_force", smallest_brute_force),
              ("long_witness", long_witness_check), ("run_extension", run_extension_check)]
    for name, function in checks:
        print(name + ": " + json.dumps(function(), sort_keys=True), flush=True)
    results = controls()
    for result in results:
        print("VERIFIED-COMPUTATIONALLY: rejected control: " + result["control"]
              + "; " + result["reason"], flush=True)
    tampered = json.loads((ROOT / "run_extension_witness.json").read_text())
    tampered["cover"].pop()
    result = expect_failure("tamper with 63-word run-extension cover",
                            lambda: run_extension_check(tampered))
    print("VERIFIED-COMPUTATIONALLY: rejected additional control: " + result["control"]
          + "; " + result["reason"], flush=True)
    print("VERIFIED-COMPUTATIONALLY: PASS; 14 exact certificates; "
          + str(sum(r["proof_nodes"] for r in checked))
          + " lower-proof nodes; 7 primary controls and 1 additional control rejected; "
            "integer arithmetic only; one process; no files written.", flush=True)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, IndexError, TypeError) as exc:
        print("REFUTED: checker rejected claim: " + str(exc), file=sys.stderr)
        sys.exit(1)
