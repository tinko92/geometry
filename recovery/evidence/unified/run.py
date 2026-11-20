"""Build in parallel, gate on public check rows, then time serially on one CPU."""
import argparse
import csv
import hashlib
import io
import json
import os
import platform
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from summarize import LABELS, measurement, schedule, summarize, validate_checks

HERE = Path(__file__).resolve().parent
BASE = "f6b03306a6b40a578428057fb42af6ea782da19d"
BASE_TREE = "29d935f3e2bf2df2e7d0103a0117187ba2d9a255191ff78f710cf755b70323b1"
TREE_PAIRS = {
    BASE_TREE: "0b4532d03b86eeb3ca5d93e4f679625008c5cdadba7ad8ac9be6659bd6b94185",
    "7863ad2022bc7e04aebe70caa806412bfff88ff2cf3ef1f1887ba82b12ea31d2":
        "b3fa66bb4590218e7e02c523d0839eaca260dfde277a4089f3887443ca16151c",
}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def tree(root):
    h = hashlib.sha256()
    for p in sorted(f for f in root.rglob("*") if f.is_file() and not f.name.startswith(".")):
        h.update(p.relative_to(root).as_posix().encode() + b"\0" + bytes.fromhex(sha(p)))
    return h.hexdigest()


def policy(label):
    if label == "old-deferred":
        return ["-DBG_BENCH_BEFORE_POLICY=deferred_nearest_output"]
    if label == "old-eager":
        return ["-DBG_BENCH_BEFORE_POLICY=eager_nearest_output"]
    if label == "new-nearest":
        return ["-DBG_BENCH_NEAREST"]
    return ["-DBG_BENCH_LEAN", "-DBG_BENCH_ULPS=" + label.split("lean")[1]]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--before", type=Path, default=Path("before"), help="old Geometry source tree")
    p.add_argument("--after", type=Path, default=Path("after"), help="new Geometry source tree")
    p.add_argument("--boost-root", type=Path, default=os.environ.get("BOOST_ROOT", "."))
    p.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    p.add_argument("--cpu", type=int, help="one allowed Linux logical CPU (required for timing)")
    p.add_argument("--seconds", type=float, default=.14)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--out", type=Path, required=True, help="new output directory")
    p.add_argument("--local-check-only", action="store_true", help="native -O0 checks; never time these binaries")
    a = p.parse_args()
    if not (.13 <= a.seconds <= .15 and a.jobs > 0):
        p.error("use 0.13..0.15 seconds, positive build parallelism")
    cxx = shutil.which(a.cxx)
    if cxx is None:
        p.error("compiler not found: " + a.cxx)
    boost = a.boost_root.resolve()
    roots = {"old": a.before.resolve() / "include", "new": a.after.resolve() / "include"}
    if not (boost / "boost/config.hpp").is_file() or not all((inc / "boost/geometry.hpp").is_file() for inc in roots.values()):
        p.error("expected Geometry source trees and umbrella BOOST_ROOT")
    affinity = []
    if not a.local_check_only:
        if (platform.system() != "Linux" or platform.machine() not in ("x86_64", "amd64")
                or a.cpu not in os.sched_getaffinity(0) or not shutil.which("taskset")):
            p.error("timing requires x86-64 Linux, taskset and one --cpu in the allowed affinity set")
        affinity = ["taskset", "-c", str(a.cpu)]
    hashes = {k: tree(inc) for k, inc in roots.items()}
    if TREE_PAIRS.get(hashes["old"]) != hashes["new"]:
        p.error(f"expected the measured historical pair or its verified Cartesian replay; got {hashes}")
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    srcdir = out / "source"
    srcdir.mkdir()
    for f, digest in {
        "simple_cartesian_benchmark.cpp": "5538591794ab26c4b93f7c58fe4f5b1aff2508ed56b802566212c8dbc0e655e1",
        "projected_benchmark_fixture.hpp": "26f2b25d6be65947825715afbe0e81a89fe46c175ed9398fa0df34fa3e044644",
    }.items():
        original = a.after / "test" / "strategies" / f
        if sha(original) != digest:
            p.error(f"{f} differs from the measured cbe774cd9 source")
        shutil.copyfile(original, srcdir / f)
    source = srcdir / "simple_cartesian_benchmark.cpp"
    flags = ["-std=c++17", "-DNDEBUG", "-fno-fast-math", "-ffp-contract=off", "-fno-omit-frame-pointer",
             "-pthread", "-DBOOST_GEOMETRY_DISABLE_DEPRECATED_14_WARNING", "-DBG_BENCH_EXACT",
             "-DBG_BENCH_COMPARE_BOUNDED"]
    flags += ["-O0"] if a.local_check_only else ["-O3", "-march=x86-64-v2"]
    args = {v: [cxx, *flags, *policy(v), "-I" + str(roots[v.split("-")[0]]), "-I" + str(boost),
                str(source), "-o", str(out / v)] for v in LABELS}
    metadata = dict(baseline=BASE, local_check_only=a.local_check_only, cpu=a.cpu, seconds=a.seconds,
                    repeats=5, affinity=affinity, platform=platform.platform(), boost_root=str(boost),
                    boost_version_sha256=sha(boost/"boost/version.hpp"),
                    compiler=subprocess.check_output([cxx, "--version"], text=True), compile_argv=args,
                    includes={k: str(inc) for k, inc in roots.items()}, tree_sha256=hashes,
                    source_sha256={f.name: sha(f) for f in srcdir.iterdir()},
                    runner_sha256=sha(Path(__file__)), summarizer_sha256=sha(HERE/"summarize.py"))
    (out / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    # Fail rather than label another copy of the old nearest implementation "new".
    guard = """#include <boost/geometry/strategies/cartesian/exact/output_policy.hpp>
#include <type_traits>
namespace r = boost::geometry::strategies::relate;
static_assert(std::is_same<r::nearest_output,r::eager_nearest_output>::value, "expected unified nearest");
"""
    subprocess.run([cxx, *flags, "-I"+str(roots["new"]), "-I"+str(boost), "-fsyntax-only", "-x", "c++", "-"],
                   input=guard, text=True, check=True)

    def build(v):
        print("Building", v, flush=True)
        with (out / (v + ".build.log")).open("w") as log:
            subprocess.run(args[v], stdout=log, stderr=subprocess.STDOUT, check=True)

    def verify_trees():
        if any(tree(inc) != hashes[k] for k, inc in roots.items()):
            raise ValueError("Geometry headers changed during campaign; use frozen exports")

    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        list(pool.map(build, LABELS))
    verify_trees()
    with (out / "replay.log").open("w") as log:
        for v in LABELS:
            cmd = [*affinity, str(out / v), "--label", v]
            subprocess.run([*cmd, "--digest-self-test", "1"], stdout=log, stderr=log, check=True)
            with (out / (v + ".check.csv")).open("w") as f:
                subprocess.run([*cmd, "--check", "1"], stdout=f, stderr=log, check=True)
        refs = validate_checks(out)
        print(f"PASS: {len(refs)} preflight cells.", flush=True)
        if a.local_check_only:
            return
        with (out / "samples.csv").open("x", newline="") as f:
            writer = None
            for v, w, op, r in schedule():
                cmd = [*affinity, str(out/v), "--label", v, "--workload", w, "--operation", op,
                       "--repeat", "1", "--repeat-offset", str(r), "--seconds", str(a.seconds)]
                data = subprocess.check_output(cmd, stderr=log, text=True)
                row, = csv.DictReader(io.StringIO(data))
                measurement(row, (v, w, op, r), refs, a.seconds)
                if writer is None:
                    writer = csv.DictWriter(f, list(row))
                    writer.writeheader()
                writer.writerow(row)
                f.flush()
    verify_trees()
    summarize(out)


if __name__ == "__main__":
    main()
