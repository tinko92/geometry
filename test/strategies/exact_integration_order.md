# Review order for the exact strategies

This is a local review series, not a proposal to change every default at
once. The base is `32cc22de5` (develop `082b99181` plus PR 1517).
The external randomized grid driver is not part of the diff. Branches
under `review/exact-final-1517/` mark the following cumulative endpoints.
Review from the preceding endpoint rather than diffing every stage against develop.
The [CPU qualification](exact_1517_cpu_results.md) summarizes the retained
implementation. The [optimization audit](exact_optimization_audit.md)
also lists rejected experiments and shortcuts for which no proof was found.

| Endpoint | Review unit |
|---|---|
| `01-native-buffer` | Ordinary buffer correctness, one cause per commit: starts at shared turns, point/line sides, represented Cartesian floating touches, skinny rings, native winding, and enabled regressions. These do not require the new exact kernel. |
| `02-strategy-hooks` | Strategy-owned turn identity, order, constructions, clipping, and output. Retain legacy fraction identity and singleton semantics for old strategies. Verify on C++14 before discussing exact dispatch. |
| `03a-arithmetic` | Bounded fixed arithmetic and scalar/point forward-error filters. Review capacities and rounding assumptions against the rational oracle. |
| `03b-construction` | Output rounding primitives, implicit intersections and operation-owned precision state. Keep this separate from algorithm changes. |
| `03c-borrowed` | Borrow original supports and compare points by certified position, not by stored segment ratios. |
| `03d-topology` | Typed predicates for bounds, generated rings, nesting and clipping. |
| `03-exact-cartesian` | Public opt-in Cartesian calls, aliasing and source lifetime; the macro-selected default pilot is a separate commit range at its end. Do not flip the global default in this PR. |
| `04-exact-optimizations` | Same semantic contract. Sterbenz dispatch, certified search limits, ring and symmetric difference shortcuts, and output copying. Keep proofs with the corresponding change. |
| `05-output-policies` | Bounded N-step construction and selective retries. Separate API review of stored, bounded and eager nearest policies; share certificate-based predicate filters. |
| `06-simple` | Uncertified Cartesian baseline, independently opt-in. |
| `07-decimal` | Authoritative scalar and nearest-even binary64 view, filtered original-source predicates. The default binary64 contract does not infer decimal intent. |
| `08-spherical` | First exact native rays and minor arcs; then directed regions, mixed-dimensional set operations and DE-9IM; finally sufficient fast paths and lossless angular export. Geographic geodesics are not included. Split public result-model review from numerical review. |
| `09-certified-buffer` | First radical signs and noding; then point/segment/line, polygon, multi-input and convex public dispatch. Exact straight offsets, arc approximation and native export have separate obligations. Refusals are transactional. |
| `09b-explicit-mp` | Orthogonal opt-in repair for existing multiprecision Cartesian buffer. Do not import a multiprecision strategy from the ordinary default headers. |
| `10-evidence` | Drivers, frozen reports and optimization audit. These data are separable from library PRs. |

Before proposing a global default, account for rounded-output collapse,
changed component counts on existing tests, named predicate and index
dispatch, mixed model support, compile cost and the C++14 default build.
For example, the one-ULP widening of a square is unequal under exact DE-9IM
but equal under the legacy polygon-equality optimization. A binary64 output
ring may also collapse even when the internal topology is exact.

For spherical overlay use the directed region to retain complement/full-sphere
and lower-dimensional results. Supplying an angular polygon destination
cannot promise lossless export of a general constructed ray. Certified
buffer results describe chordal sets with explicit budgets; this is not
a proof of universal smooth-buffer nesting for arbitrary inputs or radii.
