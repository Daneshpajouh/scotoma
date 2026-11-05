#!/usr/local/bin/python3 -I
"""Independent all-position oracle, mixed-failure controls and real CLI contracts."""
import itertools
import csv
import json
from pathlib import Path
import runpy
import unittest

HERE = Path(__file__).resolve().parent
REF = runpy.run_path(str(HERE / "reference.py"))
SUITE = runpy.run_path(str(HERE / "suite.py"))
SUPPORT = runpy.run_path(str(HERE / "run_support.py"))
run, PYTHON = SUPPORT["run"], SUPPORT["PYTHON"]


class ConformanceTests(unittest.TestCase):
    def test_exhaustive_independent_oracle(self):
        checked = 0
        for length in range(1, 7):
            for letters in itertools.product("ACGT", repeat=length):
                g = "".join(letters)
                expected = {g[:i] + g[i + 1:] for i in range(len(g))}
                self.assertEqual(REF["required_words"](g), expected)
                self.assertEqual(REF["run_count"](g), sum(1 for _ in itertools.groupby(g)))
                self.assertEqual(len(expected), REF["run_count"](g))
                self.assertEqual(SUITE["classify"](g, expected)["status"], "exact")
                checked += 1
        print(json.dumps({"exhaustive_guides": checked, "maximum_exhaustive_length": length}), flush=True)

    def test_duplicates_cannot_repair_shortfall(self):
        g = "ACGTACGT"
        words = list(REF["required_words"](g) - {g[:-1]})
        got = SUITE["classify"](g, words * 3)
        self.assertEqual(got["status"], "incomplete")
        self.assertEqual(got["missing_words"], [g[:-1]])
        self.assertGreater(got["emitted_rows"], got["r"])
        self.assertGreater(got["redundant_rows"], 0)

    def test_equal_count_invalid_substitution(self):
        g = "ACGTACGT"
        words = list(REF["required_words"](g) - {g[:-1]}) + ["N" * 7]
        got = SUITE["classify"](g, words)
        self.assertEqual(got["emitted_rows"], got["r"])
        self.assertEqual(got["status"], "invalid")
        self.assertEqual(got["missing_words"], [g[:-1]])
        self.assertEqual(got["invalid_words"], ["N" * 7])

    def test_complete_redundant(self):
        g = "AACCGGTT"
        words = sorted(REF["required_words"](g))
        got = SUITE["classify"](g, words + words[:1])
        self.assertTrue(got["word_set_exact"])
        self.assertEqual(got["status"], "redundant")
        self.assertEqual(got["redundant_rows"], 1)

    def test_terminal_run_positive_control(self):
        g = "ACGTT"
        words = {g[:i] + g[i + 1:] for i in range(len(g) - 1)}
        self.assertEqual(SUITE["classify"](g, words)["status"], "exact")

    def test_missing_entire_guide(self):
        g = "ACGT"
        got = SUITE["panel"]({"a": g, "b": g}, {"a": list(REF["required_words"](g))}, 0)
        self.assertEqual(got["guides_incomplete"], 1)
        self.assertEqual(got["per_guide"]["b"]["valid_distinct"], 0)

    def test_shared_words_keep_ownership(self):
        got = SUITE["panel"]({"a": "AAC", "b": "ACC"}, {"a": ["AC", "AA"], "b": ["CC"]}, 0)
        self.assertEqual(got["per_guide"]["b"]["missing_words"], ["AC"])

    def test_foreign_concrete_word_is_invalid(self):
        self.assertEqual(SUITE["classify"]("AAAA", ["CCC"])["status"], "invalid")

    def test_duplicate_json_key_is_rejected(self):
        """The README promises duplicate guide IDs are rejected; json.loads keeps the last
        silently, so without the hook an invalid row hides behind a valid one and the suite
        reports exact. Both directions are pinned: the duplicate is refused and the same content
        without it is still accepted, so the check discriminates rather than refusing everything.
        """
        import json as _json
        _load_json = SUITE["_load_json"]
        clean = '{"x": ["CGT", "AGT"]}'
        dup = '{"x": ["NNN"], "x": ["CGT", "AGT"]}'
        self.assertEqual(_load_json(clean), {"x": ["CGT", "AGT"]})
        self.assertEqual(_json.loads(dup), {"x": ["CGT", "AGT"]})
        with self.assertRaises(ValueError) as caught:
            _load_json(dup)
        self.assertIn("duplicate key", str(caught.exception))

    def test_invalid_inputs(self):
        for g in ("", "acgt", "ACGN", "A-CG", 4, None):
            with self.assertRaises(ValueError):
                REF["required"](g)
        for emitted in ("ACG", [None], [3]):
            with self.assertRaises(ValueError):
                SUITE["classify"]("ACGT", emitted)
        with self.assertRaises(ValueError):
            SUITE["panel"]({"a": "ACGT"}, {"unknown": []})
        with self.assertRaises(ValueError):
            SUITE["panel"]({}, {})

    def test_budget_boundaries(self):
        g = "ACGT" * 5
        self.assertTrue(REF["certificate"](g, 17)["sound_uniqueness_applies"])
        self.assertFalse(REF["certificate"](g, 18)["sound_uniqueness_applies"])
        self.assertTrue(REF["certificate"](g, 5)["paper_a_unrestricted_floor_applies"])
        self.assertFalse(REF["certificate"](g, 6)["paper_a_unrestricted_floor_applies"])
        for m in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                REF["certificate"](g, m)

    def test_cli_all_verdicts_and_malformed(self):
        directory = HERE / "build/cli_controls"
        directory.mkdir(parents=True, exist_ok=True)
        g = "ACGTACGT"
        words = sorted(REF["required_words"](g))
        for label, emitted, status, code in (
                ("exact", words, "exact", 0),
                ("incomplete", [w for w in words if w != g[:-1]], "incomplete", 1),
                ("redundant", words + words[:1], "redundant", 1),
                ("invalid", words + ["N" * 7], "invalid", 1),
                ("empty", [], "incomplete", 1)):
            path = directory / f"{label}.json"
            path.write_text(json.dumps(emitted))
            result = run(f"cli-{label}", [PYTHON, "-I", HERE / "suite.py", "--guide", g,
                                         "--emitted", path], expected=code)
            self.assertEqual(json.loads(result.stdout)["status"], status)
        path = directory / "malformed.json"
        path.write_text('{"unknown": []}')
        run("cli-unknown-id", [PYTHON, "-I", HERE / "suite.py", "--guides", HERE / "inputs/guides.tsv",
                              "--emitted", path], expected=2)
        path.write_text('{')
        run("cli-broken-json", [PYTHON, "-I", HERE / "suite.py", "--guide", g, "--emitted", path], expected=2)
        run("cli-negative-budget", [PYTHON, "-I", HERE / "reference.py", g, "--budget", "-1"], expected=2)
        count = run("cli-count-only", [PYTHON, "-I", HERE / "reference.py", g, "--count-only"])
        self.assertEqual(json.loads(count.stdout), len(words))

    def test_tsv_and_empty_word(self):
        directory = HERE / "build/cli_controls"
        directory.mkdir(parents=True, exist_ok=True)
        guides = directory / "guides.tsv"
        emitted = directory / "emitted.tsv"
        guides.write_text("guide_id\tguide\na\tA\n")
        emitted.write_text("guide_id\tpattern\na\t\n")
        result = run("cli-empty-word", [PYTHON, "-I", HERE / "suite.py", "--guides", guides,
                                      "--emitted", emitted, "--budget", "0"])
        self.assertEqual(json.loads(result.stdout)["status"], "exact")
        guides.write_text("guide_id\tguide\na\tA\na\tC\n")
        run("cli-duplicate-id", [PYTHON, "-I", HERE / "suite.py", "--guides", guides,
                                "--emitted", emitted], expected=2)

    def test_raw_alignment_adapter_rejects_corruption(self):
        adapter = runpy.run_path(str(HERE / "scoreboard.py"))["alignment_words"]
        with (HERE / "inputs/crispritz_2.9.0.tsv").open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames
            original = next(reader)
        directory = HERE / "build/cli_controls"
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "corrupted_alignment.tsv"
        for kind in ("rotated_target", "second_gap", "wrong_gap", "unknown_guide", "wrong_budget"):
            row = dict(original)
            if kind == "rotated_target":
                row["DNA"] = row["DNA"][:-3].translate(str.maketrans("ACGT", "CGTA")) + row["DNA"][-3:]
            elif kind == "second_gap":
                row["DNA"] += "-"
            elif kind == "wrong_gap":
                row["DNA"] = "-" + row["DNA"].replace("-", "")
            elif kind == "unknown_guide":
                row["crRNA"] = "A" * 20 + "NNN"
            else:
                row["Mismatches"] = "1"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
                writer.writeheader()
                writer.writerow(row)
            with self.assertRaises(ValueError, msg=kind):
                adapter(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
