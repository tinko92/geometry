"""Synthetic parser/schedule tests only; no benchmark performance measurements."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from summarize import CELLS, FIELDS, GROUPS, LABELS, schedule, summarize, validate_checks


class SummaryTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.out = Path(tmp.name)
        for v in LABELS:
            self.csv(self.out / (v + ".check.csv"), (self.row(v, w, op, 0) for w, op in CELLS))

    @staticmethod
    def row(v, w, op, r):
        return dict(label=v, workload=w, operation=op, coordinate="double", corpus_size=3,
                    calls=0, repeat=r, seconds=0, ns_per_call=0,
                    **{k: 2**63+17 if "digest" in k or "fingerprint" in k else 3 for k in FIELDS})

    @staticmethod
    def csv(path, rows):
        it = iter(rows)
        first = next(it)
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, list(first))
            writer.writeheader()
            writer.writerow(first)
            writer.writerows(it)

    def test_order_and_summary(self):
        seq = list(schedule())
        self.assertEqual(len(seq), 2925)
        self.assertEqual(len(set(seq)), 2925)
        for i in range(0, len(seq), len(LABELS)):
            window = [v for v, w, op, r in seq[i:i+len(LABELS)]]
            for g in GROUPS:
                indices = sorted(window.index(v) for v in g)
                self.assertEqual(indices[-1]-indices[0]+1, len(g))
        (self.out / "manifest.json").write_text(json.dumps(dict(seconds=.14, local_check_only=False)))
        rows = []
        for v, w, op, r in seq:
            row = self.row(v, w, op, r)
            # Repeat-dependent latency, exact 2:1 paired ratio.
            sec = (.15 + .01*r) * (2 if v.startswith("new") else 1)
            row.update(calls=3000, seconds=sec, ns_per_call=sec*1e9/3000)
            rows.append(row)
        path = self.out / "samples.csv"
        self.csv(path, rows)
        summarize(self.out)
        with (self.out / "summary.csv").open() as f:
            for row in csv.DictReader(f):
                self.assertAlmostEqual(float(row["geomean_median_paired_ratio"]),
                                       2 if row["candidate"].startswith("new") else 1)
        rows[0], rows[1] = rows[1], rows[0]
        self.csv(path, rows)
        with self.assertRaisesRegex(ValueError, "order or key"):
            summarize(self.out)

    def test_controls(self):
        path = self.out / "old-eager.check.csv"
        for name in ("input_fingerprint", "geometry_digest"):
            rows = [self.row("old-eager", w, op, 0) for w, op in CELLS]
            rows[0][name] += 1  # A single low bit would be lost in a double.
            self.csv(path, rows)
            with self.assertRaisesRegex(ValueError, "mismatch"):
                validate_checks(self.out)
        self.csv(path, rows[1:])
        with self.assertRaisesRegex(ValueError, "incomplete"):
            validate_checks(self.out)


if __name__ == "__main__":
    unittest.main()
