"""Prepare both overlays once, then run a check-only public replay per compiler."""
import argparse
import csv
import hashlib
import itertools
import json
import os
import platform
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from instrument import benchmark, predicates

HERE = Path(__file__).resolve().parent
EXACT = Path("boost/geometry/strategies/cartesian/exact")
BASE = {
    "borrowed_predicates.hpp": "1444215d06ae5a379d6c4f346b508317c4c518e0f6d5cd857417337160f257a0",
    "ulp_filter.hpp": "3596d63a3c06adadff4c4a9ce503919303dfb395320e57bc5d018704687d391a",
}
NS = (1, 2, 4, 8, 100)
VERSIONS = ("old", "new")
CELLS = set(itertools.product(
    "degenerate general decimal projected projected-local".split(),
    "union intersection difference sym_difference relate_aa relate_la relate_ll covered_by clip".split()))
KINDS = "orientation direction signed_measure unsigned_measure mixed_measure coordinate_order coordinate_equal".split()
CHECK = "coordinate corpus_digest geometry_digest corpus_size output_components output_holes output_points input_fingerprint".split()
METRICS = ("first_calls first_misses refinements refined_misses shortcut_after_miss exact_fallbacks "
           "small_integer_exact post_refine_exact tagged_exact axis_first_calls axis_first_misses "
           "axis_refined_calls axis_refined_misses").split()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_file(path, data):
    path.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")


def tree_hash(root):
    result = hashlib.sha256()
    for p in sorted(root.rglob("*")):
        if p.is_file():
            result.update((str(p.relative_to(root)) + "\0" + digest(p) + "\n").encode())
    return result.hexdigest()


def prepare(o):
    out = o.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    fixture = o.benchmark.with_name("projected_benchmark_fixture.hpp")
    meta = {"baseline": "8117188149771894eb9161521e9a8298fe3e0e98", "benchmark": str(o.benchmark.resolve()),
            "benchmark_sha256": digest(o.benchmark), "fixture_sha256": digest(fixture),
            "diagnostic_sources": {p: digest(HERE/p) for p in ("producer_diag.hpp", "instrument.py", "run.py")}}
    driver = out / "driver.cpp"
    driver.write_text(benchmark(o.benchmark.read_text()))
    shutil.copyfile(fixture, out / fixture.name)
    meta["driver_sha256"] = digest(driver)
    for v, src in zip(VERSIONS, (o.old_include, o.new_include)):
        inc = out / v / "include"
        shutil.copytree(src / "boost/geometry", inc / "boost/geometry")
        shutil.copyfile(src / "boost/geometry.hpp", inc / "boost/geometry.hpp")
        inputs = {h: digest(inc / EXACT / h) for h in BASE}
        if v == "old" and inputs != BASE:
            raise ValueError("old predicate/ULP headers are not 811718814")
        if v == "new" and any(inputs[h] == BASE[h] for h in BASE):
            raise ValueError("expected both new ULP kernels and new borrowed predicate call sites")
        original = tree_hash(inc)
        p = inc / EXACT / "borrowed_predicates.hpp"
        p.write_text(predicates(p.read_text()))
        shutil.copyfile(HERE / "producer_diag.hpp", inc / "producer_diag.hpp")
        meta[v] = {"source": str(src.resolve()), "original_headers": inputs, "original_tree_sha256": original,
                   "prepared_tree_sha256": tree_hash(inc)}
    json_file(out / "prepared.json", meta)
    print("Prepared:", out, flush=True)


def keyed(path, fields, expected):
    with path.open(newline="") as f:
        rows = {}
        for row in csv.DictReader(f):
            k = tuple(row[c] for c in fields)
            if k in rows:
                raise ValueError(f"{path}: duplicate {k}")
            rows[k] = row
    if rows.keys() != expected:
        raise ValueError(f"{path}: wrong cells or predicates: {rows.keys() ^ expected}")
    return rows


def write_csv(path, header, rows):
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def vec(values):
    return [int(values[m]) for m in METRICS]


def summarize(out, reference):
    expected = {(w, op, k) for w, op in CELLS for k in KINDS}
    fields = ["workload", "operation"]
    counts, checks, fingerprints = {}, {}, {}
    raw, pairs, totals = [], [], []
    for n, v in itertools.product(NS, VERSIONS):
        base = out / f"{v}-{n}"
        c = counts[v, n] = keyed(base / "counts.csv", fields + ["predicate"], expected)
        q = checks[v, n] = keyed(base / "check.csv", fields, CELLS)
        ref = None if reference is None else keyed(
            reference / f"{'before' if v == 'old' else 'after'}-{n}.check", fields, CELLS)
        for cell, row in q.items():
            if row["label"] != f"{v}-{n}" or row["coordinate"] != "double" or row["corpus_size"] != "3":
                raise ValueError(f"{base}: wrong corpus/label at {cell}")
            if any(float(row[m]) != 0 for m in ("calls", "seconds", "cpu_seconds", "ns_per_call")):
                raise ValueError(f"{base}: non-check replay at {cell}")
            if fingerprints.setdefault(cell[0], row["input_fingerprint"]) != row["input_fingerprint"]:
                raise ValueError(f"{base}: different inputs at {cell}")
            if ref is not None and any(row[f] != ref[cell][f] for f in CHECK):
                raise ValueError(f"{base}: uninstrumented control mismatch at {cell}")
            if v == "new" and any(row[f] != checks["old", n][cell][f] for f in CHECK):
                raise ValueError(f"N={n}: old/new public checks differ at {cell}")
        if not any(int(r["first_calls"]) for r in c.values()):
            raise ValueError(f"{base}: all probes inactive")
        raw.extend([v, n, *k, *vec(c[k])] for k in sorted(c))
    for n in NS:
        for k in sorted(expected):
            a, b = (vec(counts[v, n][k]) for v in VERSIONS)
            pairs.append([n, *k, *a, *b, *(y-x for x, y in zip(a, b))])
        for kind in KINDS:
            a, b = ([sum(int(counts[v, n][w, op, kind][m]) for w, op in CELLS)
                     for m in METRICS] for v in VERSIONS)
            totals.append([n, kind, *a, *b, *(y-x for x, y in zip(a, b))])
    cols = [v + "_" + m for v in ("old", "new", "delta") for m in METRICS]
    write_csv(out / "counts.csv", ["version", "N", *fields, "predicate", *METRICS], raw)
    write_csv(out / "paired.csv", ["N", *fields, "predicate", *cols], pairs)
    write_csv(out / "totals.csv", ["N", "predicate", *cols], totals)
    print(f"PASS: 450 three-input walks; {len(pairs)} paired predicate rows; {out/'totals.csv'}", flush=True)


def run(o):
    prep = o.out.resolve()
    data = json.loads((prep / "prepared.json").read_text())
    for v in VERSIONS:
        if tree_hash(prep / v / "include") != data[v]["prepared_tree_sha256"]:
            raise ValueError(f"{v} prepared tree changed")
    if (digest(prep / "driver.cpp") != data["driver_sha256"]
            or digest(prep / "projected_benchmark_fixture.hpp") != data["fixture_sha256"]):
        raise ValueError("prepared benchmark/fixture changed")
    boost = o.boost_root.resolve()
    if not (boost / "boost/config.hpp").is_file():
        raise ValueError("--boost-root must contain boost/config.hpp")
    cxx = shutil.which(o.cxx)
    if cxx is None:
        raise ValueError(f"compiler not found: {o.cxx}")
    out = o.results.resolve()
    out.mkdir(parents=True, exist_ok=False)
    flags = ["-std=c++17", "-O" + o.opt, "-DNDEBUG", "-fno-fast-math", "-ffp-contract=off", "-pthread",
             "-DBOOST_GEOMETRY_EXACT_FILTER_STATISTICS", "-DBG_BENCH_COMPARE_BOUNDED", "-DBG_BENCH_EXACT",
             "-DBG_BENCH_LEAN", "-DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING"]
    if o.arch:
        flags.append("-march=" + o.arch)
    json_file(out / "run.json", {"prepared": data, "platform": platform.platform(), "compiler": cxx,
              "version": subprocess.check_output([cxx, "--version"], text=True), "flags": flags,
              "boost_root": str(boost), "boost_version_header_sha256": digest(boost/"boost/version.hpp")})

    def job(job_id):
        n, v = job_id
        label = f"{v}-{n}"
        d = out / label
        d.mkdir()
        args = [cxx, *flags, f"-DBG_BENCH_ULPS={n}", "-I"+str(prep/v/"include"), "-I"+str(boost),
                str(prep/"driver.cpp"), "-o", "replay"]
        check = ["./replay", "--label", label, "--check", "1"]
        json_file(d / "commands.json", [args, check])
        print("Building", label, flush=True)
        with (d / "build.log").open("w") as log:
            subprocess.run(args, cwd=d, stdout=log, stderr=subprocess.STDOUT, check=True)
        with (d / "check.csv").open("w") as csv_out, (d / "replay.log").open("w") as log:
            subprocess.run(check, cwd=d, stdout=csv_out, stderr=log, check=True)
        print("Walked", label, flush=True)
    with ThreadPoolExecutor(max_workers=o.jobs) as pool:
        list(pool.map(job, itertools.product(NS, VERSIONS)))
    summarize(out, o.reference_dir)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    subs = cli.add_subparsers(dest="command", required=True)
    p = subs.add_parser("prepare")
    for flag in ("old-include", "new-include", "benchmark", "out"):
        p.add_argument("--" + flag, type=Path, required=True)
    r = subs.add_parser("run")
    r.add_argument("--out", type=Path, required=True, help="already prepared include trees and driver")
    r.add_argument("--results", type=Path, required=True, help="new output directory")
    r.add_argument("--boost-root", type=Path, default=os.environ.get("BOOST_ROOT", "."))
    r.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    r.add_argument("--arch", help="optional -march target (e.g. x86-64-v2)")
    r.add_argument("--opt", choices=["0", "1", "2", "3"], default="2")
    r.add_argument("--jobs", type=int, default=1)
    r.add_argument("--reference-dir", type=Path, help="also verify before-N.check / after-N.check")
    o = cli.parse_args()
    (prepare if o.command == "prepare" else run)(o)


if __name__ == "__main__":
    main()
