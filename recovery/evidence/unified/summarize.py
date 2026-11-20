"""Validate public benchmark controls and summarize matched repeat ratios."""
import argparse
import csv
import itertools
import math
from pathlib import Path
from statistics import median

NS = (1, 2, 4, 8, 100)
WORK = "degenerate general decimal projected projected-local".split()
OPS = "union intersection difference sym_difference relate_aa relate_la relate_ll covered_by clip".split()
CELLS = list(itertools.product(WORK, OPS))
TRIO = ("old-deferred", "old-eager", "new-nearest")
GROUPS = [TRIO] + [(f"old-lean{n}", f"new-lean{n}") for n in NS]
LABELS = [v for g in GROUPS for v in g]
# Each tuple is (numerator, denominator).
PAIRS = [(TRIO[2], TRIO[0]), (TRIO[2], TRIO[1]), (TRIO[1], TRIO[0])]
PAIRS += [(new, old) for old, new in GROUPS[1:]]
FIELDS = "input_fingerprint corpus_digest geometry_digest output_components output_holes output_points".split()
REPS = 5


def require(ok, message):
    if not ok:
        raise ValueError(message)


def schedule():
    for r in range(REPS):
        for i, (w, op) in enumerate(CELLS):
            for j in range(len(GROUPS)):
                k = (i + r + j) % len(GROUPS)
                g = GROUPS[k]
                if k == 0:
                    offset = (i + r) % 3
                    g = g[offset:] + g[:offset]
                    if ((i + r) // 3) % 2:
                        g = g[::-1]
                elif (i + r + k) % 2:
                    g = g[::-1]
                for v in g:
                    yield v, w, op, r


def controls(row):
    require(row["coordinate"] == "double" and row["corpus_size"] == "3", "wrong corpus/type")
    return tuple(int(row[f]) for f in FIELDS)  # Never cast uint64 fingerprints/digests to float.


def validate_checks(directory):
    refs, fingerprints = {}, {}
    for v in LABELS:
        with (directory / (v + ".check.csv")).open(newline="") as f:
            for row in csv.DictReader(f):
                w, op = row["workload"], row["operation"]
                key = v, w, op
                require(row["label"] == v and (w, op) in CELLS and key not in refs, f"bad check: {key}")
                require(row["calls"] == "0" and row["seconds"] == "0", f"not check-only: {key}")
                refs[key] = c = controls(row)
                require(fingerprints.setdefault(w, c[0]) == c[0], f"input fingerprint mismatch: {key}")
    require(len(refs) == len(CELLS)*len(LABELS), "incomplete check matrix")
    for a, b in PAIRS:
        for w, op in CELLS:
            x, y = refs[a, w, op], refs[b, w, op]
            for name, left, right in zip(FIELDS, x, y):
                require(left == right, f"{w}/{op}: {name} mismatch: {a}={left}, {b}={right}")
    return refs


def measurement(row, expected, refs, seconds):
    v, w, op, r = expected
    require((row["label"], row["workload"], row["operation"], int(row["repeat"])) == expected,
            f"wrong timing order or key; expected {expected}")
    require(controls(row) == refs[v, w, op], f"timed controls differ: {expected}")
    calls, ns, wall = int(row["calls"]), float(row["ns_per_call"]), float(row["seconds"])
    require(calls > 0 and calls % 3 == 0 and math.isfinite(wall) and wall >= seconds
            and ns > 0 and math.isclose(ns, wall*1e9/calls, rel_tol=1e-9), f"invalid window: {expected}")
    return ns


def write(path, keys, rows):
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(keys)
        writer.writerows(rows)


def summarize(directory):
    import json
    manifest = json.loads((directory / "manifest.json").read_text())
    seconds = float(manifest["seconds"])
    require(not manifest["local_check_only"] and .13 <= seconds <= .15, "not a timed server campaign")
    refs = validate_checks(directory)
    times = {}
    with (directory / "samples.csv").open(newline="") as f:
        for row, expected in itertools.zip_longest(csv.DictReader(f), schedule()):
            require(row is not None and expected is not None, "wrong number of timing rows")
            times[expected] = measurement(row, expected, refs, seconds)
    cells, summaries = [], []
    for a, b in PAIRS:
        results = []
        for w, op in CELLS:
            x, y = ([times[v, w, op, r] for r in range(REPS)] for v in (a, b))
            ratios = [left/right for left, right in zip(x, y)]
            results.append((w, op, median(ratios)))
            cells.append((a, b, w, op, median(x), median(y), median(ratios), *ratios))
        for group, key in ([("all", "all")] + [("workload", w) for w in WORK] + [("operation", op) for op in OPS]):
            selected = [ratio for w, op, ratio in results if group == "all"
                        or (group == "workload" and w == key) or (group == "operation" and op == key)]
            ratio = math.exp(sum(map(math.log, selected))/len(selected))
            summaries.append((a, b, group, key, len(selected), ratio))
            if group == "all":
                print(f"{a}/{b}: {ratio:.4f}x time ({len(selected)} cells)")
    write(directory / "paired-cells.csv", ["candidate", "reference", "workload", "operation",
          "candidate_median_ns", "reference_median_ns", "median_paired_ratio",
          *[f"ratio_rep{r}" for r in range(REPS)]], cells)
    write(directory / "summary.csv", ["candidate", "reference", "group", "key", "cells",
                                      "geomean_median_paired_ratio"], summaries)
    print(f"PASS: {len(refs)} check cells, {len(times)} timed windows. Ratios <1 are faster.")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("directory", type=Path)
    cli.add_argument("--checks-only", action="store_true")
    args = cli.parse_args()
    if args.checks_only:
        print("PASS:", len(validate_checks(args.directory)), "check cells")
    else:
        summarize(args.directory)
