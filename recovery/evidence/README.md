# Cartesian benchmark evidence (20 November 2026)

A supplement to Geometry `021b24818`, not a new benchmark or a snapshot of the library. All tables are locally generated measurements, counts or classifications. The common miniature benchmark, its projected fixture and the narrative reports already live in that revision's `test/strategies`; they are not duplicated. The labels below refer to different historical builds. Do not pool them or interpret pre-change policy names as current APIs.

| Evidence | Revisions or controls | What is retained |
|---|---|---|
| Earlier Cartesian policies | `49a77482d` vs develop `082b99181` | Three-window per-cell medians, agreement flag and ratio only for matching observables; independent synthetic XOR oracle. The historical commit may require separate access. |
| Lean N-ULP | Implementation `8c729e7c8`, driver `7169b4c2e` | Stored vs lean for N=1,2,4,8,100; 450 per-cell rows; untimed axis acceptance on original proper crossings. |
| Predicate bounds and determinant retry | `369510c05` vs `bace68a89`; retry switches in `21133c8ec` | Separate before/after bound measurements and retry off/on cells. At N=8, retry reduces source crossing points requiring refinement from 251 to 47 (614 proper crossings). These are not dynamic fallback counts. |
| Producer certificate | `811718814` vs qualified `cdfa62e2c` | Final, not preliminary, 450 public timing cells, predicate-level counts, and a deliberately difficult analytic sweep. Counts were collected separately from timing. |
| Unified certificate and eager nearest | `f6b03306a` vs `cbe774cd9` | 720 same-policy/reference comparisons, including old eager vs old deferred. Whole-operation medians and median paired ratios; full-precision replay driver. |
| Compact strategy matrix | `fab7597e9` vs `082b99181` | All 1,458 Cartesian binary64, decimal and planar buffer rows; includes unsupported, noncomparable and output-mismatch cells. |
| Mapped Cartesian workloads | GIS `021b24818`, buffer `fab7597e9`, develop `082b99181` | All 672 GIS rows; 28 native buffer pairs with a separate certified success/refusal column. File-backed drivers, reproducible feature selection and 26 operand hashes; no map coordinates. |
| Correctness context | Same-day strategy-matrix audit | Grouped verdict counts for Cartesian, Cartesian dispatch and planar buffer only. Units are policy/case ledger rows, **not** unique bugs or proof that issues were fixed. |

**Reading the tables.** Smaller time ratios favor the numerator. Keep compiler, operation, input, output policy and coordinate domain fixed. The compact suites have five workloads and nine public operations with three input pairs per cell. Timings are C++17, optimized, strict FP (no fast math or FP contraction); the primary historical toolchains were GCC 13.3 and Clang 22.1.8. Pin each paired comparison to the same available physical core. Tiny differences can be noise.

- Lean, retry, magnitude bounds and producer tables use ratios of per-cell median wall times; unified, matrix and mapped GIS use medians of same-repeat ratios. They need not equal ratios of displayed medians. The earlier policies have three windows, most later public timings five, mapped buffers six. Native absolute timings are not portable performance promises.
- The retry and lean tables keep the original develop agreement flags. They do not give develop speed ratios for unlike outputs. The retry decision can legitimately emit a different coordinate within the same ULP contract. Native controls include cardinalities, topology or spatial digests, not bitwise output equality except where explicitly stated.
- The unified historical lossless replay compared all 135 workload/input/operation transcripts across old deferred, old eager and canonical nearest on each compiler. A spatial digest match in a timing table is weaker. Exact-decimal rows use decimal literals prepared without intermediate binary64 rounding and rational output; they do not use double baselines.
- The compact matrix has 30 control mismatches with ratios withheld. Of the 672 larger GIS rows, 486 are matched candidate rows, 54 are baselines, 84 fail input validity and 48 were not scheduled. The 5,500-position Norway input was not timed. Large GIS constructive checks use component/edge digests at the recorded output grid, not an independent exact-set oracle.
- A buffer refusal is not a fast successful buffer. Certified rational output and 16-point rounded native output are different contracts and have no speed ratio. A mapped trial is planar in its stated local metric chart, not a geodesic experiment.

See [replay instructions](REPLAY.md), [selection and exclusions](SELECTION.md), and the [tables](tables). Run `python3 -B verify.py` on an unmodified copy to check file integrity and table gates. Neither publication nor verification needs network access.

**Merge-refactor control:** the six archived full-precision transcripts (two compilers, three nearest variants) have identical SHA-256 `6cc417c72373685e1ee770b2119c18fa6f6cb1423b09120e00b4b3e09b5dc9ba`. The new local replay reproduced that hash. It covers this fixed corpus, not all inputs.
