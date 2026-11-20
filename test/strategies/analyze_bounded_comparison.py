"""Validate the opt-in bounded benchmark; ratios below one mean faster."""

import csv
import math
import sys
from pathlib import Path
from statistics import median


BOUNDS = (1, 2, 4, 8, 100)
VARIANTS = ["develop"] + [f"{p}{n}" for n in BOUNDS for p in ("ulp", "lean")]
RETRY_BOUNDS = (8, 100)
RETRY_VARIANTS = ["develop"] + [f"{p}{n}-{state}" for n in RETRY_BOUNDS
                                for state in ("off", "on") for p in ("ulp", "lean")]
WORKLOADS = "degenerate general decimal projected projected-local".split()
OPS = "union intersection difference sym_difference relate_aa relate_la relate_ll covered_by clip".split()
FIELDS = "input_fingerprint corpus_digest geometry_digest output_components output_holes output_points".split()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def controls(row):
    # Keep digests as integers, never floating point (some exceed 2^53).
    require(row["coordinate"] == "double" and row["corpus_size"] == "3",
            "unexpected coordinate type or corpus size")
    return tuple(int(row[f]) for f in FIELDS)


def check(directory, variants):
    refs, fingerprints = {}, {}
    for variant in variants:
        with (directory / f"{variant}-check.csv").open() as f:
            for row in csv.DictReader(f):
                w, op = row["workload"], row["operation"]
                key = (variant, w, op)
                require(row["label"] == variant and w in WORKLOADS and op in OPS
                        and key not in refs, f"unexpected check row {key}")
                refs[key] = values = controls(row)
                require(fingerprints.setdefault(w, values[0]) == values[0],
                        f"native input fingerprint mismatch: {key}")
                peers = []
                if variant.startswith("lean"):
                    peers.append("ulp" + variant[4:])
                if variant.endswith("-on"):
                    peers.append(variant[:-3] + "-off")
                for peer in peers:
                    for field, a, b in zip(FIELDS, refs[peer, w, op], values):
                        require(a == b, f"{key}: {peer} {field} mismatch: {a} != {b}")
    require(len(refs) == len(variants) * len(WORKLOADS) * len(OPS), "incomplete check matrix")
    with (directory / "input-fingerprints.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["workload", "input_fingerprint"])
        writer.writerows((w, fingerprints[w]) for w in WORKLOADS)
    return refs


def timed(directory, refs, repetitions, workloads, variants):
    samples = {}
    with (directory / "results.csv").open() as f:
        for row in csv.DictReader(f):
            key = (row["label"], row["workload"], row["operation"])
            require(key in refs and key[1] in workloads, f"unexpected timed row {key}")
            require(controls(row) == refs[key], f"timed/check fingerprint mismatch: {key}")
            r, calls, ns = int(row["repeat"]), int(row["calls"]), float(row["ns_per_call"])
            group = samples.setdefault(key, {})
            require(0 <= r < repetitions and r not in group, f"invalid repeat {r}: {key}")
            require(calls > 0 and calls % 3 == 0 and math.isfinite(ns) and ns > 0,
                    f"invalid measurement: {key}, repeat {r}")
            group[r] = ns
    for v in variants:
        for w in workloads:
            for op in OPS:
                require(len(samples.get((v, w, op), {})) == repetitions,
                        f"incomplete timing: {(v, w, op)}")
    with (directory / "comparison.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        if variants == RETRY_VARIANTS:
            writer.writerow(["N", "layout", "workload", "operation", "develop_ns", "off_ns", "on_ns",
                             "on_over_off", "off_over_develop", "on_over_develop", "develop_output_matches_off"])
            for n in RETRY_BOUNDS:
                for p, layout in (("ulp", "stored"), ("lean", "lean")):
                    for w in workloads:
                        for op in OPS:
                            keys = [(v, w, op) for v in ("develop", f"{p}{n}-off", f"{p}{n}-on")]
                            d, off, on = (median(samples[k].values()) for k in keys)
                            writer.writerow([n, layout, w, op, d, off, on, on/off, off/d, on/d,
                                             refs[keys[0]] == refs[keys[1]]])
            return
        writer.writerow(["N", "workload", "operation", "develop_ns", "stored_ns", "lean_ns",
                         "stored_over_develop", "lean_over_develop", "lean_over_stored",
                         "develop_output_matches_stored"])
        for n in BOUNDS:
            for w in workloads:
                for op in OPS:
                    keys = [(v, w, op) for v in ("develop", f"ulp{n}", f"lean{n}")]
                    d, s, l = (median(samples[k].values()) for k in keys)
                    writer.writerow([n, w, op, d, s, l, s/d, l/d, l/s, refs[keys[0]] == refs[keys[1]]])


if __name__ == "__main__":
    try:
        directory = Path(sys.argv[1])
        args = sys.argv[2:]
        variants = RETRY_VARIANTS if args[:1] == ["--det-retry"] else VARIANTS
        if variants == RETRY_VARIANTS:
            args = args[1:]
        refs = check(directory, variants)
        if args[0] == "timed":
            timed(directory, refs, int(args[1]), args[2:], variants)
        else:
            require(args == ["check"], "expected check or timed REPETITIONS WORKLOAD...")
    except (ValueError, OSError, KeyError, IndexError) as exc:
        sys.exit(f"Bounded comparison failed: {exc}")
