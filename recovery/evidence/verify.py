"""Check the published package and optional locally regenerated GIS/buffer inputs."""
import argparse
import csv
import hashlib
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INDEX = "SHA256SUMS"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(name):
    with (ROOT / "tables" / (name + ".tsv")).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def verify(fixtures=None):
    expected = {}
    for line in (ROOT / INDEX).read_text().splitlines():
        digest, name = line.split("  ", 1)
        assert name not in expected
        expected[name] = digest
    paths = list(ROOT.rglob("*"))
    assert not any(p.is_symlink() for p in paths)
    actual = {p.relative_to(ROOT).as_posix(): sha(p) for p in paths if p.is_file() and p != ROOT / INDEX}
    assert expected == actual, "unexpected, missing, or modified file; use an unmodified export"

    counts = {"earlier-policies":720, "lean":450, "retry":360, "producer":450, "unified":720,
              "construction-counts":5, "retry-counts":10, "public-counts":70, "crossing-counts":50,
              "crossing-sweep":10, "predicate-bounds":270, "matrix":1458, "gis":672,
              "mapped-buffer":28, "fixture-hashes":26, "verdict-counts":314}
    assert all(len(rows(name)) == n for name,n in counts.items())
    for name, metric in (("matrix", "median_paired_time_ratio"), ("gis", "paired_time_ratio")):
        for r in rows(name):
            eligible = r["status"] in ("baseline", "matched")
            assert bool(r[metric]) == eligible
            if eligible: assert 0 < float(r[metric]) < math.inf
    assert Counter(r["status"] for r in rows("gis")) == dict(matched=486, baseline=54, **{"invalid-input":84, "not-scheduled":48})
    assert sum(r["status"]=="output-mismatch" for r in rows("matrix"))==30
    for r in rows("mapped-buffer"):
        assert (r["certified_status"]=="success") == bool(r["certified_success_us"])
        assert (r["certified_status"]=="cannot_certify") == bool(r["certified_refusal_us"])
    for r in rows("earlier-policies"):
        assert bool(r["ratio_if_matched"]) == (r["observable_matches"]=="True")
    if fixtures:
        gis, buffers = fixtures
        for r in rows("fixture-hashes"):
            if r["group"]=="gis":
                p = gis/"inputs"/r["sample"]/(r["operand"]+".wkt")
            else:
                p = buffers/"fixtures"/(r["sample"]+".wkt")
            assert sha(p) == r["sha256"], r["sample"]
    print(f"PASS: {len(expected)} indexed files, {sum(counts.values())} table rows"
          + (", 26 regenerated operand hashes" if fixtures else ""))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fixtures", nargs=2, type=Path, metavar=("GIS_OUTPUT","BUFFER_OUTPUT"))
    verify(p.parse_args().fixtures)
