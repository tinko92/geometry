# Exact versus simple Cartesian: overhead profile

Measurement only: no strategy or algorithm implementation was changed.
The corpus and benchmark call body are unchanged from
`simple_cartesian_performance.md`.

## Method

CPU-only on the same Xeon Gold 6530 server, GCC 13.3, CPU 2. Native timings use
the existing optimized binaries, three serial rotated repetitions per group,
0.1 seconds per repetition. Separate Callgrind 3.22 runs use the same optimization
and floating-point flags plus debug information: C++17, `-O3 -g -DNDEBUG
-ffp-contract=off -fno-omit-frame-pointer -pthread`.

The profiling wrapper includes the original benchmark source and calls its
unchanged `execute` function. Parsing, input validation, warm-up and final harness
teardown are outside instrumentation; each operation's setup and destruction
remain inside its measured call. Each profile measures 500 complete three-input
batches: 1,500 calls. There are 54 profiles and 81,000 measured calls.

The server prohibits unprivileged perf and the account lacks sudo rights.
Installed gprofng timer pilots were unreliable and were discarded. Valgrind was
downloaded directly on the server and extracted in the experiment directory;
there was no system package installation or kernel setting change.

**Callgrind measures executed instructions, not cycles or proportional elapsed
time.** In particular, the attribution percentages below are not wall-time
percentages. They identify added computational work; they do not establish a
recoverable speedup or a native cache/memory-latency breakdown.

Validation: every exclusive instruction sum matches both the raw totals footer
and the standard annotation tool. All recorded execute-call counts are 1,500.
All 54 fixed-call checksums match native execution, and all 54 complete output
controls under instruction emulation match native controls. Sources and headers
are fingerprinted. No correctness-relevant changes were made for profiling.

## Relative costs

The primary comparison excludes the two previously identified invalid simple
symmetric-difference groups. Those profiles are retained separately. For the
25 matching groups, the geometric mean exact/simple ratio is **1.974 in native
time** and **2.125 in executed instructions**. The full 27-group native ratio is
2.069, but includes the two non-equivalent simple outputs.

| Operation | Groups | Native Time Ratio | Instruction Ratio |
| --- | ---: | ---: | ---: |
| Union | 3 | 1.889 | 2.006 |
| Intersection | 3 | 1.680 | 1.776 |
| Difference | 3 | 1.984 | 2.033 |
| Symmetric difference, degenerate category only | 1 | 1.679 | 1.719 |
| Relate areal/areal | 3 | 1.990 | 2.156 |
| Relate linear/areal | 3 | 2.364 | 2.598 |
| Relate linear/linear | 3 | 2.115 | 2.428 |
| Point covered-by | 3 | 1.820 | 1.865 |
| Line/polygon intersection | 3 | 2.138 | 2.428 |

## Where the extra work goes

For an equal-call-count mixture of the 25 matching groups, mean instructions per
call are 64,050 exact versus 31,615 simple. This pooled ratio differs from the
geometric mean above. The following is a disjoint, source-based attribution of
the net extra instructions, using source-line positions and function owners,
not occurrences of exact types in template arguments:

| Added Work | Share of Extra Instructions |
| --- | ---: |
| Predicate/construction paths and associated numerical support | 43.5% |
| Point/support representation, including variant dispatch | 29.1% |
| Allocation/deallocation bookkeeping | 10.3% |
| Fixed-capacity arithmetic and backend checks | 4.4% |
| Ring topology/parenting | 2.3% |
| Remaining algorithms, standard library, explicit copying and unresolved code | 10.4% |

The main files are `forward_error.hpp`, `borrowed_predicates.hpp`,
`borrowed_bounds.hpp`, `borrowed_records.hpp`, and variant dispatch. Allocation
calls average 72.4 malloc calls per operation for exact versus 53.8 simple.
These include algorithmic containers and metadata, not just numerical caches.
Nested allocator layers are not added together when counting allocations.

Some useful contrasts:

- Decimal linear/linear relate spends almost the same instructions in the
  boundary-checker body in both configurations. Additional work is concentrated
  in endpoint lookup, lexicographic comparison, Root visitation and scalar
  filtering, rather than simply more traversal.
- Decimal point covered-by executes about 461 instructions/call exact versus
  237 simple. There are no allocator or fixed-integer routines in that profile:
  native winding, certified orientation and dispatch already expose overhead.
- Fixed arithmetic and its checked-backend validation together account for only
  2.2% of exact instructions across the primary mixture. They are not the general
  explanation for the approximately twofold slowdown.
- Decimal difference is an exception: backend capacity/limb checks alone cost
  9,438 instructions/call, 7.33% of that case's exact work. The backends are
  statically sized, such as 257/515/1032-bit integers, not dynamically sized
  `cpp_int` allocations.
- BGL-source instructions contribute about 0.015% of the net added work. This
  does not measure the entire shared graph-processing cost, including library
  callbacks and allocation, and is not a BG-versus-GEOS comparison.

Inline attribution and unresolved libc routines limit the semantic precision
of the buckets. An unidentified allocator helper is attributed through its
exclusively allocator incoming call graph; other unresolved routines remain
separate. Disassembly identifies copy/clear routines among that remainder.

## Best next targets

1. **Known native/search-coordinate comparisons.** Finite binary64 search
   metadata currently goes through the general scalar-filter machinery even
   though it has zero error. Specializing that case can preserve exact order
   while retaining the current infinity, integral and implicit-point paths.
   See `borrowed_bounds.hpp:330` and `forward_error.hpp:561`.
2. **Repeated operand decoding.** Share transient coordinate extraction and
   visitation across the X/Y parts of lexicographic comparison and endpoint
   lookup. This need not change the borrowed representation or retain new point
   copies. See `borrowed_predicates.hpp:295`.
3. **Allocation and parenting scratch.** Inspect operation-local reservations
   and repeated boundary/residual construction. Do not infer from the flat
   allocator total that all of it belongs to parenting or numerical caches.
4. **Forced floating-point spills (subsequently removed).** At profiling time,
   `arithmetic.hpp:84-87` used volatile
   stores/loads for every basic rounded operation. The current supported
   contraction/reassociation/precision contract may make some redundant. These
   lines cost about 926 instructions per primary-case call on average, including
   the necessary arithmetic itself, so this is not an isolated removable-cost
   estimate. The qualified removal and new timings, including nearest-output
   modes, are recorded in `exact_no_volatile_performance.md`.
5. **Checked fixed-backend operations.** Review redundant internal checks only
   where existing capacity proofs cover every intermediate. Keep checked input
   admission and all required exact fallbacks.

Replacing robust orientation with a raw determinant is not an optimization
candidate: it changes the required guarantees. The other prospective speedups
above have not been measured by ablation.

## Reproducing

`profile_simple_cartesian_benchmark.sh` builds symbol-rich fixed-call drivers,
collects instructions, and exports standard Callgrind function tables. Supply
`PROFILE_ROOT`, `BOOST_ROOT`, `EXACT_INCLUDE` and `VALGRIND_ROOT`; optional
`PROFILE_ROUNDS` defaults to 500 and `TASKSET_CPU` to 2. Run its `build`, `collect`
and `reports` actions in sequence. Analyze without launching any target:

```sh
python3 test/strategies/analyze_simple_cartesian_profile.py DATA_DIRECTORY \
  --medians native-analysis.json --output analysis.json --text analysis.txt
```

Raw working results: `/tmp/bg-compact-profile-20261116`. Server experiment:
`/home/tinko/dev/bg-exact-compact-20261116/profile`. Compressed profiles,
validation, fresh timings and analysis are retained under the umbrella checkout's
`.geometry-plans/exact-compact-profile-20261116`. No commits or pushes.
