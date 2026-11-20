# Cartesian experiment recovery

This is the Cartesian work from the 2026-11-20 experiment, not a proposed
Boost.Geometry patch series. Start at `d50cc9053`: history through native
overlay/relate, decimal, simple and the earlier output policies is retained
with its original commit IDs. The subsequent spherical implementation is
omitted. All 44 commits after the spherical section up to `021b24818` were
replayed; see [original and backup IDs](cartesian-commits.tsv) when a historical
report cites the mixed branch. Some historical reports mention both domains.

All non-spherical production headers at this backup tip are byte-identical
to the source experiment. The two strategies test manifests and the standalone
qualification script omit the absent experimental spherical tests. Ordinary
upstream spherical code and tests are unchanged. No build products, server
access configuration or external randomized-grid checkout are needed to
recover the Cartesian source.

The [portable evidence supplement](evidence/README.md) includes additional
per-operation timings, diagnostics and reproduction recipes. It does not
redistribute external datasets or issue attachments. Verify its manifest
with `python3 -B recovery/evidence/verify.py`.

For a standalone smoke test from the Geometry root, set BOOST_ROOT to a local
Boost checkout with headers and use a fresh output directory:

```sh
CASE_FILTER=lean_bounded_ulp sh test/strategies/run_exact_primitives.sh /tmp/bg-cartesian-smoke
```
