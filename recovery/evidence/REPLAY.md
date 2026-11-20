# Replay from local sources

Use temporary output directories, an assembled Boost header tree (`BOOST_ROOT`), a C++17 compiler (`CXX`), and a Geometry repository containing the named Git objects (`GEOMETRY`). Set `PACK` to this evidence folder, `WORK` to an empty scratch directory. The snippets use shell variables, not fixed installation locations. Git archives read objects without changing a checkout.

On the public backup this supplement lives at `recovery/evidence`. Historical IDs throughout these notes are the **old** IDs. Resolve them with [the adjacent replay map](../cartesian-commits.tsv) (old SHA, replayed SHA, message). The backup preserves the 07-decimal ancestry and replays 44 later commits. Prefer the replayed ID when mapped, otherwise require the original object to exist. Work directory names below remain the old short IDs for clarity. Do not substitute the tip for a missing experimental object.

```sh
resolve() {
    local mapped
    mapped=$(awk -F '\t' -v old="$1" 'index($1,old)==1 {print $2}' "$PACK/../cartesian-commits.tsv") || return 1
    git -C "$GEOMETRY" rev-parse --verify "${mapped:-$1}^{commit}"
}
for rev in 021b24818 7169b4c2e 21133c8ec 811718814 cdfa62e2c f6b03306a cbe774cd9 fab7597e9 082b99181; do
    resolved=$(resolve "$rev") || exit 1
    mkdir "$WORK/$rev"
    git -C "$GEOMETRY" archive "$resolved" include test/strategies | tar -x -C "$WORK/$rev"
done
cd "$WORK"
```

**Source-crossing diagnostics.** Compile `drivers/construction_counts.cpp` with both include and test directories from `7169b4c2e`; compile `retry_counts.cpp` the same way with `21133c8ec`. Run without arguments; compare CSV values to the corresponding TSV. Neither program times anything. For example:

```sh
"$CXX" -std=c++17 -O3 -DNDEBUG -fno-fast-math -ffp-contract=off \
  -I7169b4c2e/include -I7169b4c2e/test/strategies -I"$BOOST_ROOT" \
  "$PACK/drivers/construction_counts.cpp" -o construction
./construction
```

For the timing experiment use the tracked `run_simple_cartesian_benchmark.sh` from the corresponding driver revision with `BENCH_COMPARE_BOUNDED=1` or `BENCH_COMPARE_DET_RETRY=1`. Follow that version's header configuration and timing switches. The magnitude experiment uses the `369510c05` driver and a 0.12-second window, independently of retry. The November 19 table is historical evidence, not a rerun using the current nearest policy.

**Producer-bound experiment.** The public runner instruments disposable include copies, not installation headers. A predicate-level refinement is a retry decision, not an allocation or an individual coordinate operation. Nested predicates can be counted separately. For targeted counts only positive first-call rows are retained; the two compilers' counted rows agreed.

```sh
python3 -B "$PACK/producer/run.py" prepare --old-include 811718814/include \
  --new-include cdfa62e2c/include --benchmark 811718814/test/strategies/simple_cartesian_benchmark.cpp --out probe
python3 -B "$PACK/producer/run.py" run --out probe --results counts \
  --cxx "$CXX" --boost-root "$BOOST_ROOT" --opt 3
```

Compile `drivers/crossing_sweep.cpp` uninstrumented with each version's include root; run with a budget of 1, 2, 4, 8 or 100. Every pass checks four predicates for each of 1,537 analytic near-collinear crossings. Each measurement runs for at least 0.15 seconds, reporting ns per sign; seven alternating measurements per arm. For counter output define `BOOST_GEOMETRY_EXACT_FILTER_STATISTICS` and use the corresponding instrumented include tree and the companion `producer_diag.hpp`.

To reproduce the uninstrumented **public** producer experiment on Linux x86-64, name the two extracted trees `before` and `after` in scratch (or make local symbolic links to the extracts). Run `bash "$PACK/producer/time.sh" "$CXX" repeat "$CPU"`, then `python3 -B "$PACK/producer/summarize.py" repeat`. The summarizer checks complete five-repeat coverage and equality of the input, output and cardinality controls before publishing ratios. These times are not mixed with the targeted sweep.

**Unified nearest.** The archived driver at `cbe774cd9` differs from the one at `021b24818` in policy-selection preprocessing. The supplied runner reads that exact version from the after tree and checks the driver and fixture hashes. Select one available CPU and use a new output directory for each compiler. Timed mode requires Linux x86-64 and taskset; `--local-check-only` is an untimed build/validation mode. The summarizer's synthetic tests run with `python3 -B "$PACK/unified/test_summary.py"`.

```sh
python3 -B "$PACK/unified/run.py" --before f6b03306a --after cbe774cd9 \
  --boost-root "$BOOST_ROOT" --cxx "$CXX" --cpu "$CPU" --jobs 2 --out unified-run
```

For a full-precision transcript compile `drivers/full_precision.cpp` with `-DBG_BENCH_EXACT`, the common `cbe774cd9/test/strategies` driver on the include path, and the tested version's headers. For the old build use `-DBG_BENCH_BEFORE_POLICY=deferred_nearest_output` or `=eager_nearest_output`; for the new build use `-DBG_BENCH_NEAREST -DEXPECT_CANONICAL_NEAREST`. Compare the 135 output lines byte for byte. The independent `xor_oracle.cpp` uses seeded triangles and sheared boxes; compile against `021b24818` and run with no arguments. It compares exact vertical and horizontal occupancy profiles and boundary multiplicity.

The matrix's independent Python oracle is available as `drivers/slab_oracle.py`; its self-test uses synthetic squares only. Its `compare(op, a, b, result)` consumes valid, closed polygonal WKT in exact rational arithmetic, tests open-slab membership of regularized sets and declines large inputs. It neither validates arbitrary WKT nor tests lower-dimensional output. The issue corpus is deliberately not part of the package.

**Strategy matrix.** The binary64 driver is the tracked `simple_cartesian_benchmark.cpp` at `fab7597e9`. The extra drivers here are Cartesian-only extracts of the actual decimal and compact buffer drivers. Create the decimal fixture header locally:

```sh
python3 -B "$PACK/recipes/decimal_inputs.py" fab7597e9/test/strategies/simple_cartesian_benchmark.cpp decimal-fixtures.hpp
"$CXX" -std=c++17 -O3 -DNDEBUG -fno-fast-math -ffp-contract=off -DDECIMAL -DCAND \
  -Ifab7597e9/include -I. -I"$BOOST_ROOT" "$PACK/drivers/decimal.cpp" -o decimal
./decimal decimal intersection 0
```

For the rational baseline omit `-DCAND` and use develop headers. The compact buffer driver uses `point`, `bent-line` or `parcel` as its first argument, then `buffer` and a repeat index; `-DCAND` selects certified output. Build the unmodified native branch against each tree for the native comparison. Rotate variants for five windows of at least 0.13 s; compare inputs and full observed outputs before calculating ratios.

**File-backed cases.** Use only input files you are entitled to access. Recreate GIS inputs with `python3 -B "$PACK/recipes/gis.py" "$GEOMETRY" gis-data --revision "$(resolve 021b24818)"` (reads immutable Git objects). For the buffer recipe export the two source files from the pinned revision:

```sh
mkdir map-source
git -C "$GEOMETRY" archive "$(resolve fab7597e9)" example/data/roads.wkt test/algorithms/buffer/data/gr_ll.wkt | tar -x -C map-source
python3 -B "$PACK/recipes/buffer.py" map-source buffer-data
python3 -B "$PACK/verify.py" --fixtures gis-data buffer-data
```

Compile `drivers/gis.cpp` with `-DEXACT -DNEAREST` against `021b24818` for new nearest, no policy flags against develop, or `-DRATIONAL -DFILTERED` vs `-DRATIONAL` for decimal. Arguments: operand directory, `polygon` or `track`, operation, `check` or `time`, grid exponent (16 for country windows, 28 for degree-like units). Always run check before time, match controls, then rotate five windows of at least 0.15 s. Example: `./gis gis-data/inputs/it-window polygon intersection check 16`.

Native buffer builds use no macro; use `-DCERTIFIED` for the rational-output candidate. Arguments: `road` or `coast`, metric WKT file, signed offset, tolerance, batch count (zero for preflight). Use +25 m for roads, +/-250 m for coasts. The measured tolerance was the absolute offset divided by 32 (0.78125 m or 7.8125 m). Run six rotated batches, using the same batch size for the paired native calls.
