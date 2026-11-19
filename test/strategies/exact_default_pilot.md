# Opt-in Cartesian Exact Default Pilot

This pilot changes only selected public calls that **omit** a strategy.
No macro means no change. To opt in, compile the entire program with
`BOOST_GEOMETRY_ENABLE_CARTESIAN_EXACT_DEFAULT_PILOT` defined, C++17 or later,
and strict IEEE binary64 floating-point evaluation: no fast-math,
reassociation, excess precision, or implicit contraction; rounding to nearest,
ties to even; gradual underflow (no FTZ or DAZ). For example, use
`-fno-fast-math -ffp-contract=off` on GCC or Clang. Do not define the macro in
a single source file or mix different language/FP modes or macro settings
across translation units that use Boost.Geometry. This header-only feature
changes inline template definitions. The macro attests to the FP flags:
preprocessor checks alone cannot detect every incompatible flag.

C++14 and detectable fast-math, finite-only, excess-precision and
non-binary64 builds keep the legacy route without loading the exact header
through the selector. Non-nearest runtime rounding also takes the legacy
route. The pilot does not change the rounding mode or FP flags.

## Admitted Calls

Use two concrete `model::polygon<model::d2::point_xy<double>>` inputs (default
clockwise, closed, standard storage), valid as represented. No correction,
validation or snapping is done. Every ordinate, including those in holes,
must be an integer multiple of `2^-64` in `[-2^16, 2^16]`. The deliberately
small grid bounds native differences and products well inside normal
binary64 range: differences are multiples of `2^-64` of magnitude at most
`2^17`; products have grid `2^-128` and magnitude at most `2^34`.
This is an admission restriction, not a replacement for the existing
exact-kernel filter proofs and fallbacks.

| Call without strategy | Additional requirement |
| --- | --- |
| `relation(a, b)` | None |
| `relate(a, b, mask)` | `de9im::mask`, not a static mask |
| `intersection`, `union_`, `difference`, `sym_difference` | Destination is exactly `model::multi_polygon<model::polygon<model::d2::point_xy<double>>>` |

All other signatures keep the legacy path. Selection happens before dynamic
visitation, so a variant containing an otherwise admitted polygon is still
legacy. Calls with any explicit strategy, including `bg::default_strategy{}`,
are unchanged. The shared Cartesian umbrella, insert helpers, named
predicates, validation, index, unary operations, and buffer are not opted in.
For example, two polygons one ULP apart can be disjoint according to pilot
`relation` yet touching according to legacy `touches`. Use an explicit
strategy for all related calls when they must agree.

## Native Output And Feedback

The pilot is equivalent to explicitly passing
`strategies::relate::cartesian_exact<double, strategies::relate::approximate_output>{}`
on admitted calls. Internal exactness does not imply exact, nearest-even,
or topology-preserving native export. The overlay emits the current stored
approximation without output-only refinement. Distinct exact points can
collapse; for example #663 can append a zero-area ring even though the exact
input relation has overlapping interiors. A later operation sees only the
exported doubles, not the original construction records. Validate native
output before reusing it and do not expect a chain of native-output calls to
preserve exact set identities. Use explicit rational output when exact
constructed coordinates must survive a call.

## Focused Test

`exact_default_pilot.cpp` is registered in CMake and Boost.Build. Standalone
examples from the Geometry root (use separate executables):

```sh
CXX=clang++
BOOST_ROOT=/path/to/boost
SRC=test/strategies/exact_default_pilot.cpp
PILOT=-DBOOST_GEOMETRY_ENABLE_CARTESIAN_EXACT_DEFAULT_PILOT
ON=-DBOOST_GEOMETRY_TEST_EXPECT_EXACT_DEFAULT_PILOT
OFF=-DBOOST_GEOMETRY_TEST_EXPECT_LEGACY_DEFAULT_PILOT
$CXX -Iinclude -isystem "$BOOST_ROOT" -std=c++17 -O1 -ffp-contract=off $PILOT $ON "$SRC" -o /tmp/bg-pilot-on
$CXX -Iinclude -isystem "$BOOST_ROOT" -std=c++17 -O1 -ffp-contract=off $OFF "$SRC" -o /tmp/bg-pilot-off
$CXX -Iinclude -isystem "$BOOST_ROOT" -std=c++14 -O1 -ffp-contract=off $PILOT $OFF "$SRC" -o /tmp/bg-pilot-cxx14
$CXX -Iinclude -isystem "$BOOST_ROOT" -std=c++17 -O1 -ffast-math $PILOT $OFF "$SRC" -o /tmp/bg-pilot-fast
/tmp/bg-pilot-on
/tmp/bg-pilot-off
/tmp/bg-pilot-cxx14
/tmp/bg-pilot-fast
```

The active test compares all four overlays against explicit approximate
output, including pre-populated destinations that own both input polygons.
It also compares a crossing whose approximate export differs from nearest
and legacy output, and checks the excluded signatures and fallback paths.
