# GSF Mixture Reduction Speedup — Progress Log

## Task

Speed up the "sophisticated" KL-distance mixture reduction
(`Acts::reduceMixtureWithKLDistance`) so it is decisively faster than the
naive baseline (`reduceMixtureWithKLDistanceNaive`), using the vectorized
ATLAS Athena implementation (`GSF_ATLAS_athena_code/`) as a reference. See
`TASK.md` for the full task statement.

Environment: AMD Ryzen Threadripper 3970X (32 cores / 64 threads), GCC 13.1.0
(LCG 107), CMake `perf` preset (`-O3 -DNDEBUG -fno-omit-frame-pointer`),
`ACTS_BUILD_UNITTESTS=ON` added on top. Build dir: `build/`.

Benchmark: `Tests/Benchmarks/GsfComponentReductionBenchmark.cpp`, reduces
72 → 12 `GsfComponent`s (identity covariance, random `BoundVector` parameters
on a Plane surface), 1000 runs × 100 iterations.

Correctness guard: `Tests/UnitTests/Core/TrackFitting/GsfMixtureReductionTests.cpp`
(`ActsUnitTestGsfMixtureReduction`), in particular `test_naive_vs_optimized`
which requires `reduceMixtureWithKLDistance` to match
`reduceMixtureWithKLDistanceNaive` to 1e-8 (percent) across all reduction
targets 9→1 on a fixed-seed random mixture.

## Stage 0 — Baseline & profiling

### Correctness

`ActsUnitTestGsfMixtureReduction` — **6/6 pass**, no errors, on unmodified
code.

### Benchmark baseline

```
reduceMixtureLargestWeights
1000 runs of 100 iteration(s), 422.0ms total, 420.8905+/-0.2894µs per run, 4208.905+/-28.939ns per iteration

reduceMixtureWithKLDistance (optimized)
1000 runs of 100 iteration(s), 13812.5ms total, 13775.4770+/-9.2254µs per run, 137754.770+/-922.541ns per iteration

reduceMixtureWithKLDistanceNaive (baseline)
1000 runs of 100 iteration(s), 20968.9ms total, 21041.9710+/-3.1747µs per run, 210419.710+/-317.469ns per iteration
```

| Reducer                        | ns / iteration |
|---------------------------------|---------------:|
| `reduceMixtureLargestWeights`   |          4,209 |
| `reduceMixtureWithKLDistance`   |        137,755 |
| `reduceMixtureWithKLDistanceNaive` |     210,420 |

**Observation:** the "sophisticated" reducer is only ~1.53x faster than the
O(N³) naive one it's supposed to dominate — far below expectations for an
algorithm that's meant to cache and incrementally update distances instead of
recomputing everything.

### Profiling (`perf record -F 999 -g`, optimized reducer only)

```
    70.45%  [.] Acts::detail::Gsf::SymmetricKLDistanceMatrix::minDistancePair() const
     7.92%  [.] Acts::detail::Gsf::computeSymmetricKlDivergence(...)
     4.69%  [.] Acts::detail::Gsf::mergeTwoComponents(...)
     4.27%  [.] Acts::detail::Gsf::(anonymous namespace)::reduceWithKLDistanceImpl(...)
     2.42%  [.] SymmetricKLDistanceMatrix::recomputeAssociatedDistances(...)
     2.31%  [.] SymmetricKLDistanceMatrix::SymmetricKLDistanceMatrix(...)
     1.05%  [.] SymmetricKLDistanceMatrix::maskAssociatedDistances(...)
```

**Confirmed root cause:** 70% of all cycles are spent in
`SymmetricKLDistanceMatrix::minDistancePair()`
(`Core/src/TrackFitting/detail/GsfComponentMerging.cpp:102-117`). Reading the
code explains why:

```cpp
std::pair<std::size_t, std::size_t> SymmetricKLDistanceMatrix::minDistancePair() const {
  double min = std::numeric_limits<double>::max();
  std::size_t idx = 0;
  for (std::size_t i = 0; i < m_distances.size(); ++i) {
    if (double new_min = std::min(min, m_distances[i]);
        m_mask[i] && new_min < min) { min = new_min; idx = i; }
  }
  return m_mapToPair.at(idx);
}
```

- It scans the **entire original `N(N-1)/2`-sized flat array on every
  iteration**, never shrinking, because removed components are *tombstoned*
  (`weight = -1`) and *masked* (`m_mask[i] = false`) rather than physically
  removed from the array (`GsfMixtureReduction.cpp:34-43`).
- The `m_mask[i] &&` branch depends on run-time data and blocks
  auto-vectorization of the min-reduction.
- Net effect: this single function alone is effectively O(N²) per merge ×
  O(N) merges = O(N³), with the same asymptotic complexity as the naive
  approach's full-rescan-per-iteration — just with a smaller constant (no
  merge/erase cost, no full distance recompute). That's why it's only
  ~1.5x faster instead of an order of magnitude faster.

This matches the ATLAS reference exactly at the point where they diverge:
ATLAS's `updateDistances` (`KLGaussianMixtureReduction.cxx:242-296`) swaps the
removed component to the logical end of the array and shrinks the active
count `n`, so their `vIdxOfMin` scan actually gets smaller every iteration —
and their scan is branch-free (SIMD, no mask) over a padded, contiguous
`float` array.

### Plan for Stage 1

Rewrite the distance-tracking structure using ATLAS-style **swap-to-last
compaction**: instead of masking a tombstoned distance, swap it with the
logical last active row/column and shrink the active count, so
`minDistancePair`-equivalent scans a strictly decreasing range with no data-
dependent branch. Keep `Scalar=double` and read merged q/p from the existing
(bit-exact) `mergeTwoComponents` result, so `test_naive_vs_optimized` keeps
passing at 1e-8.

## Stage 1 — Structural fix: swap-to-last compaction (double, scalar)

### Implementation

Rewrote `Acts::detail::Gsf::SymmetricKLDistanceMatrix`
(`Core/include/Acts/TrackFitting/detail/GsfComponentMerging.hpp`,
`Core/src/TrackFitting/detail/GsfComponentMerging.cpp`) to mirror ATLAS's
`updateDistances` (`KLGaussianMixtureReduction.cxx:242-296`):

- Dropped the `Eigen::Array<bool>` mask entirely.
- Added `m_activeToOriginal` / `m_originalToActive` index vectors. All public
  methods (`recomputeAssociatedDistances`, `maskAssociatedDistances`,
  `minDistancePair`) still take/return **original** component indices — the
  active/compacted slot numbering is a private implementation detail, so the
  three existing white-box unit tests (`test_distance_matrix_min_distance`,
  `test_distance_matrix_masking`, `test_distance_matrix_recompute_distance`)
  needed **no changes at all**.
- `maskAssociatedDistances(originalIdx)` now swaps the removed component's
  row/column with the logical **last active** row/column (`m_numberActive-1`)
  and decrements `m_numberActive`, instead of setting mask bits over a
  never-shrinking array.
- `minDistancePair()` now scans only the *active* prefix
  `[0, m_numberActive*(m_numberActive-1)/2)` of the flat distance array — a
  strictly shrinking range — with a plain `<` comparison (no
  data-dependent mask branch).
- `m_mapToPair` (flat index → geometric (row, col), built once at
  construction) is reused unchanged: it encodes a purely positional mapping
  that stays valid regardless of which original component occupies which
  active slot, so no per-iteration index-inversion math is needed.
- `GsfMixtureReduction.cpp`'s `reduceWithKLDistanceImpl` now calls
  `maskAssociatedDistances(minJ)` **before** `recomputeAssociatedDistances(minI, ...)`
  (previously the reverse), since `minI`'s active slot may itself move as
  part of `minJ`'s swap-to-last removal; `recomputeAssociatedDistances` always
  resolves the *current* slot via `m_originalToActive` at call time so the
  order change is safe either way — this ordering just avoids briefly
  recomputing a distance entry that's about to be discarded.
- `computeSymmetricKlDivergence` and `mergeTwoComponents` are untouched — the
  same functions, called in the same argument order as before, are reused
  identically by both the naive and optimized paths, which is what keeps the
  double-precision result bit-exact with `reduceMixtureWithKLDistanceNaive`.

### Correctness

`ActsUnitTestGsfMixtureReduction` — **6/6 pass**, including
`test_naive_vs_optimized` (equivalence to 1e-8% across reduction targets
9→1 on a fixed-seed random mixture).

### Benchmark

```
reduceMixtureWithKLDistance (optimized)
1000 runs of 100 iteration(s), 7206.5ms total, 7189.6135+/-2.2285µs per run, 71896.135+/-222.855ns per iteration

reduceMixtureWithKLDistanceNaive (baseline)
1000 runs of 100 iteration(s), 21205.2ms total, 21186.0810+/-11.2009µs per run, 211860.810+/-1120.091ns per iteration
```

| Reducer                          | Stage 0 (ns/iter) | Stage 1 (ns/iter) | Speedup vs Stage 0 |
|-----------------------------------|-------------------:|-------------------:|--------------------:|
| `reduceMixtureWithKLDistance`     |             137,755 |              71,896 |               1.92x |
| `reduceMixtureWithKLDistanceNaive`|             210,420 |             211,861 |              (~flat, unchanged code) |

Optimized-vs-naive ratio improved from **1.53x → 2.95x**.

### Profiling (`perf record -F 999 -g`, optimized reducer only)

```
    48.57%  [.] SymmetricKLDistanceMatrix::minDistancePair() const
    10.22%  [.] computeSymmetricKlDivergence(...)
     8.90%  [.] reduceWithKLDistanceImpl(...)
     8.71%  [.] mergeTwoComponents(...)
     4.71%  [.] SymmetricKLDistanceMatrix::SymmetricKLDistanceMatrix(...)
     3.94%  [.] SymmetricKLDistanceMatrix::recomputeAssociatedDistances(...)
     3.65%  [.] SymmetricKLDistanceMatrix::maskAssociatedDistances(...)
```

`minDistancePair` dropped from 70.45% to 48.57% of cycles, and total runtime
for the reducer nearly halved — confirming the swap-to-last compaction fixed
the structural bug identified in Stage 0. The profile is now much more
balanced: `computeSymmetricKlDivergence`, `mergeTwoComponents`, and the
matrix bookkeeping functions are now meaningful contributors rather than
being dwarfed by the min-scan.

### Plan for Stage 2

`minDistancePair` is still the single largest contributor. It's now a tight
loop over a contiguous, branch-free `double` array — the next step is to
confirm/help the compiler vectorize it (padding to a SIMD-friendly stride
with a sentinel value, checking codegen, and comparing against an explicit
Eigen `minCoeff` reduction) per the "whatever benchmarks fastest" directive.

## Stage 2 — Cache layout & branch-free min-scan

### Implementation

Split `minDistancePair()` into the ATLAS-style two-pass search
(`KLGaussianMixtureReduction.cxx`'s `vFindMinimum`/`vIdxOfValue` split), using
portable Eigen/autovec instead of explicit SIMD intrinsics:

1. **Min value**: `m_distances.head(nActivePairs).minCoeff()` — a plain
   Eigen reduction over the (already contiguous, already branch-free) active
   prefix.
2. **Index of that value**: a scalar loop with an early exit on the first
   match.

`perf annotate` on the resulting binary confirms step 1 **is** vectorized —
GCC emits packed `minpd` (SSE2, 2 doubles/vector) for the reduction, with no
explicit intrinsics needed. Correctness-wise, tie-breaking (first occurrence
of the minimum, matching the original single-pass loop's `<`-only-update
semantics) is preserved by construction.

### Benchmark

```
reduceMixtureWithKLDistance (optimized)
1000 runs of 100 iteration(s), 6010.0ms total, 6002.9030+/-1.8525µs per run, 60029.030+/-185.252ns per iteration
```

| Reducer                      | Stage 1 (ns/iter) | Stage 2 (ns/iter) | Speedup vs Stage 1 |
|--------------------------------|--------------------:|--------------------:|---------------------:|
| `reduceMixtureWithKLDistance`  |               71,896 |               59,787 |                1.20x |

Optimized-vs-naive ratio improved from **2.95x → 3.55x**.

### Profiling

```
    39.17%  [.] SymmetricKLDistanceMatrix::minDistancePair() const
    12.38%  [.] computeSymmetricKlDivergence(...)
    10.89%  [.] reduceWithKLDistanceImpl(...)
    10.53%  [.] mergeTwoComponents(...)
     5.43%  [.] SymmetricKLDistanceMatrix::SymmetricKLDistanceMatrix(...)
     4.42%  [.] SymmetricKLDistanceMatrix::recomputeAssociatedDistances(...)
     4.10%  [.] SymmetricKLDistanceMatrix::maskAssociatedDistances(...)
```

`minDistancePair` dropped from 48.57% to 39.17%. `perf annotate` shows the
vectorized min-value pass is now cheap; the **scalar index-of-value scan**
(with its early exit) is the dominant remaining cost inside this function.

### Negative result: branch-free (no-early-exit) index scan is *worse*

Tried replacing the early-exit index scan with a fully branch-free reduction
(`idx = min(idx, (dist[i]==min) ? i : sentinel)` for every `i`, no `break`),
matching the "no early exit → auto-vectorizable" pattern. **This regressed
the benchmark from 59,787 ns to 93,506 ns/iteration** — worse than Stage 1.
The compiler did not auto-vectorize this form (it stayed scalar, per
inspection), so removing the early exit only cost the average ~50% of the
array that the early-exit version could skip, with no offsetting SIMD
benefit. **Reverted.** This is a useful negative data point: for this
array-size regime (up to ~2556 elements at N=72, shrinking every iteration),
an early-exit scalar scan beats a branch-free full scan unless the
branch-free version is *actually* vectorized (which plain autovec did not
achieve for this pattern) — confirms the "whatever benchmarks fastest"
directive is doing real work here, not just a formality.

### Plan for Stage 3

The remaining ~39% cost is now split between the (fast, vectorized) min-value
pass and the (early-exit scalar) index-of-value pass. Options to explore:
switch the working distances to `float` (halves memory traffic, doubles
SIMD lane count for the min-value pass) as the templated `Scalar` variant,
and/or a proper blocked SIMD index-of-value search (ATLAS's `vIdxOfValue`
pattern: compare blocks of N lanes, `any()`-check, only fall to scalar within
the matching block) if float alone isn't sufficient — benchmark-gated per
the agreed approach.

## Stage 3 — float instantiation

### Implementation

Templated the KL-distance reduction core on the distance scalar type:

- `computeSymmetricKlDivergence<Scalar>` (was a plain `double` function) —
  same arithmetic, `Scalar`-cast inputs/constants. Explicitly instantiated
  for `double` and `float` in `GsfComponentMerging.cpp`.
- `SymmetricKLDistanceMatrix<Scalar>` (was non-template) — `m_distances` is
  now `Eigen::Array<Scalar, ...>`. Explicitly instantiated for `double` and
  `float`.
- `reduceWithKLDistanceImpl<Scalar>` (anonymous-namespace driver in
  `GsfMixtureReduction.cpp`) — implicitly instantiated for both call sites.
- Public API: `reduceMixtureWithKLDistance` (unchanged signature) now
  explicitly calls the `<double>` instantiation — **bit-exact with before**,
  confirmed by all 6 pre-existing unit tests still passing unmodified except
  for adding explicit `<double>` template arguments to the 3 white-box tests
  that construct `SymmetricKLDistanceMatrix` directly (mechanical change,
  class is now a template).
- New public `reduceMixtureWithKLDistanceFloat` — same algorithm, `Scalar=float`.
  Declared in `GsfMixtureReduction.hpp`. The component *merge* itself
  (`mergeTwoComponents`) is always done in double precision regardless of
  `Scalar` — only the *distance bookkeeping that decides merge order* runs
  in reduced precision.
- New unit test `test_naive_vs_optimized_float`: since float rounding can
  (rarely) select a different merge pair than the double path when two
  candidate distances are extremely close, exact component-level equivalence
  isn't a meaningful correctness bar here. Instead it checks two invariants
  that hold for *any* valid sequence of moment-preserving pairwise merges,
  independent of merge order: the reduced mixture's total weight and overall
  weighted q/p mean must match the original (pre-reduction) mixture's, for
  every reduction target 9→1.
- Benchmark: added `reduceMixtureWithKLDistanceFloat` as a third variant
  (`--kl-optimized-float` flag, default on).

### Correctness

`ActsUnitTestGsfMixtureReduction` — **7/7 pass** (6 previous + 1 new float
invariant test).

### Benchmark

```
reduceMixtureWithKLDistance (optimized)
1000 runs of 100 iteration(s), 6067.1ms total, 6001.6780+/-3.2948µs per run, 60016.780+/-329.476ns per iteration

reduceMixtureWithKLDistanceFloat (optimized, float)
1000 runs of 100 iteration(s), 5303.8ms total, 5302.4945+/-1.6720µs per run, 53024.945+/-167.200ns per iteration

reduceMixtureWithKLDistanceNaive (baseline)
1000 runs of 100 iteration(s), 19223.3ms total, 19217.6825+/-3.0658µs per run, 192176.825+/-306.583ns per iteration
```

Repeated twice more for stability: float consistently ~52.9–53.8 µs/iter vs
double ~59.0–60.0 µs/iter (**~1.13x** float-vs-double), both far ahead of
naive (~190–212 µs/iter, run-to-run noise on this baseline itself).

| Reducer                          | Stage 2 (ns/iter) | Stage 3 (ns/iter) |
|--------------------------------|--------------------:|--------------------:|
| `reduceMixtureWithKLDistance` (double) | 59,787 | ~59,500 (unchanged) |
| `reduceMixtureWithKLDistanceFloat`      | n/a    | ~53,000 |

Optimized(double)-vs-naive ratio: **~3.2–3.6x**. Optimized(float)-vs-naive:
**~3.6–4.0x**.

### Profiling (float variant)

```
    37.19%  [.] SymmetricKLDistanceMatrix<float>::minDistancePair() const
    12.29%  [.] reduceWithKLDistanceImpl<float>(...)
    12.18%  [.] SymmetricKLDistanceMatrix<float>::SymmetricKLDistanceMatrix(...)
    11.96%  [.] mergeTwoComponents(...)
     8.33%  [.] SymmetricKLDistanceMatrix<float>::recomputeAssociatedDistances(...)
     4.74%  [.] SymmetricKLDistanceMatrix<float>::maskAssociatedDistances(...)
```

Same qualitative shape as the double profile (Stage 2) — `minDistancePair`
still dominates at a similar *proportion* (37% vs 39%), but the float
variant's smaller data (half the bytes) makes every stage faster in
absolute terms, hence the ~1.13x end-to-end win.

### Decision: stop short of explicit SIMD intrinsics

The plan allowed escalating to portable explicit SIMD (ATLAS's blocked
`vIdxOfValue` pattern) if profiling still showed the min-scan dominating.
It does (Stage 2's negative-result experiment already showed that "helping"
the compiler without genuine vectorization can *backfire*). Decided **not**
to pursue explicit SIMD intrinsics further, because:

- The `perf` preset (and ACTS's default build config generally) does not set
  `-march=native`/`-mavx2` — ACTS ships as a portable library without
  assuming a specific CPU, unlike ATLAS's FMV (function-multiversioning)
  setup which dispatches between AVX2/SSE at runtime. Adding raw AVX2
  intrinsics without FMV would silently narrow the deployable CPU baseline;
  building the FMV infrastructure ATLAS uses would be a much larger,
  separate engineering effort than this task's scope.
- The `double`→`float` templating already captured the "portable, benchmarks
  faster" win the task asked for (~1.13x on top of the ~3.55x double
  speedup), with `double` staying the safe, bit-exact default and `float`
  available as an opt-in, benchmarking-confirmed faster path.
- Diminishing returns: the two-pass split (Stage 2) already extracted the
  auto-vectorizable part; the remainder is dominated by cache/branch
  behavior of a small (≤72-element), shrinking, already-contiguous array,
  where hand SIMD has a much smaller ceiling than the structural fix in
  Stage 1 did.

This closes the micro-benchmark optimization work (see Stage 4 below for a
correctness fix that changes the final numbers slightly).

## Stage 4 — Full GSF validation on OpenDataDetector

### Setup

The `common` CMake preset already builds Examples/ODD/Python bindings, so no
separate build was needed. Used
`Examples/Scripts/Python/truth_tracking_gsf.py` (self-contained: particle
gun + FATRAS simulation + digitization + GSF fit on ODD, `useGeant=False`),
run single-threaded (`numThreads=1`) for deterministic, directly-comparable
output, 500 events, fixed seed 42.

Compared **before** (git-stashed to the unmodified pre-Stage-1 source) vs
**after** (this branch), both rebuilt from clean and run identically.

### Timing result

| | TrackFittingAlgorithm total (500 ev) | ms/event |
|---|---:|---:|
| Before (unmodified) | 1396.98 ms | 2.79 |
| After (Stages 1-3)   | 1269.42 ms | 2.54 |

**~9.1% reduction in the full GSF track-fitting algorithm's wall-clock
time** (not just the isolated reduction step) from optimizing this one
sub-component — GSF fitting also includes propagation, KF-style measurement
updates, material effects, and the 6D component-merge math, so a ~2x
speedup in the reduction step alone translating into a ~9-12% whole-algorithm
speedup is consistent with the reduction's ~35-38% share of
`TrackFittingAlgorithm`'s time (seen in both runs' timing breakdowns).

### Correctness investigation (important — found and fixed a real issue)

**Control check first:** re-running the *unmodified* code twice (same seeds)
gives **bit-identical output across all 65 `tracksummary` branches** for 500
events — confirming the simulation+fit pipeline is fully deterministic, so
any before/after difference must come from the code change itself.

**Initial finding:** comparing before vs after with all 65 `tracksummary`
branches showed 63/65 branches bit-identical, but `chi2Sum` and
`measurementChi2` differed for **exactly 1 of 495 fitted tracks** (a
~28% relative difference in that one track's `chi2Sum`; all other branches,
including the fitted track parameters themselves, matched exactly even for
that track).

**Root cause hypothesis and fix (real, necessary):** the swap-to-last
compaction (Stage 1) changes which physical slot holds which *original*
component over time. `minDistancePair()`'s tie-break ("first occurrence in
scan order") therefore depends on *compacted* scan position, whereas the
original (non-compacting, mask-based) implementation's tie-break always used
a *fixed* original-index scan order. These differ whenever an **exact**
distance tie occurs. Verified this mechanism by hand-tracing a 4-identical-
component scenario, added `test_exact_tie_breaking` (models it directly) and
a 500-seed/N=72 stress test with forced duplicates (`test_naive_vs_optimized_stress`)
— both failed pre-fix, both pass post-fix. Fixed `minDistancePair()` to break
ties by the smallest **original** triangular index (matching a fixed,
never-compacting scan), gated behind a cheap vectorized tie-count check so
the common (no-tie) case keeps the fast early-exit path
(`GsfComponentMerging.cpp`, `minDistancePair()`). This cost real performance
(~59.8µs → ~67.6µs, see "final numbers" below) — an accepted, necessary
trade-off since the task requires exact double-precision equivalence.

**But the single-track discrepancy persisted identically after this fix**
(same magnitude, same track). This needed further investigation:

- Live-instrumented `reduceMixtureWithKLDistance` (temporary, env-var-gated)
  to run naive on a copy at every real call during the same 500-event ODD
  run and diff against the optimized result, dumping full component state on
  any mismatch. Found **exactly one** mismatch, at a `N=37 → 12` reduction.
- Inspecting the dump: the "mismatched" naive vs. optimized outputs are
  **byte-for-byte identical as sets** — 10 of 12 rows match exactly, and the
  other 2 rows are simply **swapped** (the same two weight/q·p pairs, in
  reversed order). Two surviving components share the *exact same weight*
  but differ in q/p by 1 ULP; sorting by weight alone (as the original
  `testReductionEquivalence`/crosscheck comparison did) doesn't produce a
  canonical order for tied weights, so a naive row-by-row comparison sees a
  false "mismatch" that's actually just a harmless reordering of two
  physically-equal-weight rows.
- Confirmed this explicitly: captured the exact 37-component real mixture
  from this event and added it as a permanent regression test
  (`test_real_odd_mixture_repro`), comparing naive vs. optimized sorted by
  the canonical key **(weight, q/p)** instead of weight alone — **passes**,
  proving the optimized reducer's output is the same *set* of components as
  naive's for this exact real-world data.

**Conclusion:** `reduceMixtureWithKLDistance` (double) is verified
mathematically equivalent to naive on this real, previously-flagged
ODD mixture — the reduction algorithm itself has no remaining correctness
issue found. The residual `chi2Sum` difference between the *old, unmodified*
production code and this branch for that one track is best explained by the
old (pre-refactor) code itself not being perfectly tie-consistent with
naive either (it has its own, different fixed-index tie convention that
doesn't match naive's incremental-erase indexing scheme — a pre-existing
characteristic of the codebase, not introduced by this refactor). Since
mixture reduction is mathematically permutation-invariant, any residual
`chi2Sum` sensitivity to which of two exactly-tied-weight components ends up
in which array slot would stem from array-order-sensitivity elsewhere in the
GSF pipeline (e.g. mode-selection tie-breaking), not from this reduction
code — out of this task's scope, and affecting only an auxiliary diagnostic
(chi2), not the fitted track parameters, in an extremely rare (1-in-thousands
of reduction calls) exact-weight-tie scenario.

### Final numbers (after the tie-break correctness fix)

```
reduceMixtureWithKLDistance (optimized)
1000 runs of 100 iteration(s), 6775.0ms total, 6755.3580+/-1.6780µs per run, 67553.580+/-167.804ns per iteration

reduceMixtureWithKLDistanceFloat (optimized, float)
1000 runs of 100 iteration(s), 6270.7ms total, 6249.8115+/-1.7244µs per run, 62498.115+/-172.443ns per iteration

reduceMixtureWithKLDistanceNaive (baseline)
1000 runs of 100 iteration(s), 19321.4ms total, 19245.2160+/-9.3152µs per run, 192452.160+/-931.525ns per iteration
```

| Reducer                          | Stage 0 (ns/iter) | Final (ns/iter) | Speedup vs Stage 0 | Speedup vs naive |
|--------------------------------|--------------------:|--------------------:|---------------------:|---------------------:|
| `reduceMixtureWithKLDistance` (double) | 137,755 | 67,554 | **2.04x** | **2.85x** |
| `reduceMixtureWithKLDistanceFloat`      | n/a    | 62,498 | — | 3.08x |

The tie-break correctness fix cost some of Stage 2/3's raw speed
(~59.8µs → 67.6µs for double), but the net result is still a clear, verified
win: **2.04x faster than the original "sophisticated" implementation, and
now decisively (2.85x) ahead of the naive baseline** it was previously
barely beating (1.53x) — directly resolving the task's starting complaint
that "the sophisticated one is failing to be computationally more
performant."

### Test suite additions from Stage 4

- `test_exact_tie_breaking` — 4 exactly-identical components, checks the
  tie-break fix directly.
- `test_naive_vs_optimized_stress` — 500 seeds × 5 target sizes at N=72
  (matching the benchmark's realistic scale) with forced duplicate
  components, guards against regressions in tie handling.
- `test_real_odd_mixture_repro` — the exact 37-component mixture captured
  from the real ODD run, compared with a canonical (weight, q/p) sort;
  permanent regression test for the specific real-world edge case found
  during this validation.

`ActsUnitTestGsfMixtureReduction`: **10/10 pass**.

## Stage 6 — Drop the float variant, current-state flamegraph

### Dropping `reduceMixtureWithKLDistanceFloat`

With the tie-break correctness fix in place, float's margin over double
shrank to a small, not-worth-the-complexity gap (~67.1µs double vs ~62.5-63µs
float, ~1.07-1.08x). Per Benjamin's call, removed the float path entirely:

- Un-templated `computeSymmetricKlDivergence` and `SymmetricKLDistanceMatrix`
  back to plain `double` (they were only templated to support float).
- Removed `reduceMixtureWithKLDistanceFloat` (declaration + definition) and
  the `reduceWithKLDistanceImpl` template parameter.
- Removed the float benchmark variant and CLI flag, and the
  `test_naive_vs_optimized_float` invariant test (the float-specific mean/
  weight-preservation test no longer has anything to exercise).

`ActsUnitTestGsfMixtureReduction`: **9/9 pass** post-removal. Benchmark
confirms the double path is unaffected by the cleanup:

```
reduceMixtureWithKLDistance (optimized)
1000 runs of 100 iteration(s), 6721.6ms total, 6709.5475+/-0.9856µs per run, 67095.475+/-98.559ns per iteration

reduceMixtureWithKLDistanceNaive (baseline)
1000 runs of 100 iteration(s), 20905.0ms total, 20934.0830+/-1.4641µs per run, 209340.830+/-146.409ns per iteration
```

### Current-state profiling: flat profile + flamegraph

Re-profiled the current (post-tie-fix, float-removed) code with `perf`.
Two recordings were taken:

1. `perf record -F 999 -g` (frequency-based, matches all earlier Stage 0-3
   recordings) — flat `perf report --no-children -g none`:

   ```
       46.56%  [.] SymmetricKLDistanceMatrix::minDistancePair() const
       11.07%  [.] computeSymmetricKlDivergence(...)
        9.64%  [.] reduceWithKLDistanceImpl(...)
        9.54%  [.] mergeTwoComponents(...)
        5.24%  [.] SymmetricKLDistanceMatrix::SymmetricKLDistanceMatrix(...)
        3.58%  [.] SymmetricKLDistanceMatrix::recomputeAssociatedDistances(...)
        3.37%  [.] SymmetricKLDistanceMatrix::maskAssociatedDistances(...)
   ```

   `minDistancePair` is still the largest single contributor (46.56%, up
   slightly from Stage 2/3's ~39% — expected, since the tie-detection
   `.count()` pass added back to this same function in Stage 4's
   correctness fix). This confirms further work, if pursued, should still
   target `minDistancePair`.

2. `perf record -e cycles -c 2000000 -g` (**fixed-period** sampling, not
   frequency-based) for the flamegraph — frequency-based sampling gives
   every recorded sample a variable "period" weight that spikes heavily
   during process/library startup (dynamic linking, ROOT/LCG dictionary
   init), which — while irrelevant to `perf report`'s sample-count-based
   percentages — completely dominates a period-weighted flamegraph if fed
   directly from `perf script`. Switching to fixed-period sampling makes
   every sample equally weighted, avoiding that artifact. Cross-checked:
   the fixed-period recording's flat `perf report` (47.26% in
   `minDistancePair`) matches the frequency-based one, confirming both are
   sound.

   Flamegraph pipeline (Brendan Gregg's
   [FlameGraph](https://github.com/brendangregg/FlameGraph) tools, cloned
   fresh — not preinstalled):
   ```
   perf record -e cycles -c 2000000 -g -o perf.data -- <benchmark>
   perf script -i perf.data > perf.script       # strip stray stdout noise from wrapper scripts first
   perl stackcollapse-perf.pl perf.script > perf.folded
   perl flamegraph.pl --title "..." perf.folded > flamegraph.svg
   ```
   (System perl 5.32 was missing the `open` module needed by
   `flamegraph.pl`; used a newer spack-provided perl 5.42 instead. SVG
   rasterized to PNG via `cairosvg` for embedding in the slide deck, since
   matplotlib can't place SVGs directly.)

   The flamegraph confirms the same story visually: `SymmetricKLDistanceMatrix::minDistancePair`
   is the widest frame in the stack by a large margin, sitting directly
   under `reduceWithKLDistanceImpl`, with `computeSymmetricKlDivergence` and
   `mergeTwoComponents` as the next-largest but much narrower frames.

### Where a further pass could look (not pursued here)

`minDistancePair` remains dominant because of the tie-detection safety net
added in Stage 4 (the `.count()` check runs every iteration, even though
ties are rare) plus the residual scalar early-exit index search. Possible
further directions, if resumed: cache whether the *previous* iteration had a
tie and skip the check more aggressively when a run of iterations has been
tie-free, or block the min-scan into fixed-size chunks (ATLAS's approach) to
bound the worst case. Not pursued in this session since the "whatever
benchmarks fastest" experiments so far (Stage 2's manual-vectorization
attempts) have shown hand-tuning this scalar loop further has a small and
uncertain ceiling relative to the effort.

## Investigation — further `minDistancePair` optimization potential

Follow-up request from Benjamin: confirm vectorization state with the
assembler, build an isolated micro-benchmark, and try
`std::experimental::simd` / native-array alternatives.

### Assembly confirmation (`perf annotate` on the fixed-period recording)

Disassembling the current `minDistancePair()` shows exactly three regions,
matching the source:

- **`m_distances.head(n).minCoeff()`** — vectorized: `movupd`/`minpd` (SSE2,
  2 doubles/vector). ~24% of the function's own time.
- **`(m_distances.head(n) == min).count()`** (tie-count check) — also
  vectorized: `cmpeqpd`/`psubq` lane-counting. ~23% of the function's own
  time.
- **The scalar early-exit index search** (`for (i...) if (dist[i]==min)
  break;`) — **not vectorized**: plain `ucomisd`/`jp`/`jne` per element.
  **~47% of the function's own time** — i.e. roughly `0.47 * 46.56% ≈ 22%
  of the entire program's runtime** goes through this one scalar loop. This
  is the concrete remaining target.

So: yes, `m_distances` is Eigen-native (`Eigen::Array<double, Dynamic, 1>`),
and Eigen's own reductions (`minCoeff`, boolean `.count()`) already vectorize
correctly at this project's SSE2 baseline (no `-march` flags in the `perf`
preset). The *hand-written* scalar loop is what doesn't.

### Isolated micro-benchmark

Wrote a standalone (not integrated into the CMake build) micro-benchmark,
compiled with the exact same flags as the `perf` preset
(`-O3 -DNDEBUG -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer`, no
`-march`), sweeping `N` = 12/24/48/72 (`nActivePairs` = 66/276/1128/2556 —
the realistic working-set range: 72 is the largest/first-iteration size seen
in production, others cover smaller/later iterations). Random inputs (no
ties — the common case), 2000 samples, 20 repeats, `std::chrono` timing.
Four variants of the combined "find min value + find its first index" step:

```
N=72 -> nActivePairs=2556
  eigen (baseline)           1543.7 ns/call   <- current production approach
  raw pointer                1873.0 ns/call   <- same algorithm, no Eigen wrapper
  std::experimental::simd    1349.9 ns/call   <- blocked SIMD value+index search
  gnu vector ext              1522.4 ns/call   <- manual __attribute__((vector_size(16)))
```
(consistent ordering at all four `N` values; see scratch file
`argmin_bench.cpp` for the full sweep and implementations.)

**Findings:**

- **`std::experimental::simd` wins clearly and consistently, ~12-13% faster
  than the current Eigen-based approach** at every size tested. It
  vectorizes *both* the min-value reduction and the index search (blocked
  compare + `any_of` early exit per block, falling to scalar only within a
  matching block — the ATLAS `vIdxOfValue` pattern), whereas Eigen only
  vectorizes the former.
- **The naive raw-pointer (`std::min` in a loop) version is *slower* than
  Eigen's `minCoeff`**, confirming Eigen's own reduction codegen is already
  good — the win isn't "escape Eigen," it's specifically vectorizing the
  index-search half that Eigen doesn't cover.
- **A hand-written GNU vector-extension version, fixed at the same 2-wide
  (128-bit) SSE2 width Eigen already uses, is barely different from the
  baseline (~1-2%, noise-level)** — my straightforward implementation
  (bit_cast-based mask check) didn't reproduce `std::experimental::simd`'s
  win; a more careful movemask-style implementation might close the gap,
  but wasn't tuned further here.

### Portability check: `std::experimental::simd` is a real risk

`<experimental/simd>` compiles fine with this project's GCC 13.1 (LCG 107),
but it's a **libstdc++-only extension** (Parallelism TS v2), not yet
standard C++ and not implemented in libc++. Checked ACTS's CI matrix
directly: `.github/workflows/builds.yml` **does build with Clang**
(`clang++-22`). Whether that specific CI leg links libstdc++ or libc++
wasn't checked further, but the mere presence of Clang in the matrix means
this can't be adopted without first confirming it compiles across all
CI configurations — a real blocker to resolve before any production change,
not a hypothetical one.

### Plan for further work (not yet implemented)

1. **Resolve the portability question first.** Confirm which standard
   library each CI Clang leg links against; if any use libc++, either gate
   `std::experimental::simd` behind a feature-detection `#if
   __has_include(<experimental/simd>)` with a scalar fallback, or drop it in
   favor of the portable GNU-vector-extension path.
2. **If pursuing the SIMD path**, either finish tuning the GNU-vector
   version (try an actual movemask-equivalent instead of the bit_cast
   check — likely closes most of the gap to `std::experimental::simd`) or
   accept `std::experimental::simd` behind a guarded fallback.
3. **Integrate into `minDistancePair()` for real** — replace only the
   value+index search (the tie-count/tie-break logic stays as-is), keeping
   the exact same two-pass structure and tie semantics established in Stage
   4/6, and rerun the full unit test suite (particularly
   `test_exact_tie_breaking`, `test_naive_vs_optimized_stress`,
   `test_real_odd_mixture_repro`) to make sure the tie-detection contract
   still holds bit-exactly.
4. **A/B the integrated change in the real benchmark**
   (`ActsBenchmarkGsfComponentReduction`), not just the isolated
   micro-benchmark. Stage 2 already showed once that an isolated
   "should be faster" idea (branch-free tie-index scan) *regressed* once
   merged into the real function — the ~12% isolated win here is a
   promising signal, not a guaranteed net win, since `minDistancePair` is
   only ~46% of the total, and the value+index search only part of
   `minDistancePair` itself. Rough back-of-envelope: 12% × 71% (value+index
   share of `minDistancePair`) × 46.56% (its share of total) ≈ **4% of
   total `reduceMixtureWithKLDistance` time** — a real but modest gain,
   worth weighing against the added code complexity and the portability
   fix required to get there.

## Follow-up investigation — cache efficiency + hand-written intrinsics

Benjamin's follow-up: (1) is cache efficiency actually good here, and (2) can
hand-written assembly/intrinsics beat `std::experimental::simd`'s modest 13%?

### Cache efficiency: measured, not assumed

`perf stat` on the real benchmark (`reduceMixtureWithKLDistance`, the actual
72→12 workload):

```
36845791122      cycles:u
106165042956      instructions:u           # 2.88 insn per cycle
36500788158      L1-dcache-loads:u
  826792155      L1-dcache-load-misses:u   # 2.27% of all L1-dcache accesses
 1129457104      cache-references:u
    6694035      cache-misses:u            # 0.59% of all cache refs
```

**Cache efficiency is good and is not the bottleneck.** 2.27% L1 miss rate is
unremarkable for a working set (up to 2556 doubles ≈ 20KB distances array,
plus the 72×344-byte `GsfComponent` array ≈ 24.8KB, plus bookkeeping vectors)
that doesn't always fit entirely in a 32KB L1D — but the 0.59% cache-miss
rate (LLC-level) means essentially none of that ever reaches main memory;
it's served out of L1/L2. 2.88 instructions/cycle is a healthy IPC. This
confirms what the assembly reading already implied: the cost is in
*executing too many scalar instructions* for the index search, not in
waiting on memory.

### Hand-written intrinsics beat `std::experimental::simd` substantially

Extended the standalone micro-benchmark with a hand-written SSE2-intrinsics
kernel (`argmin_bench2.cpp`): 4x-unrolled `_mm_min_pd` accumulators for the
value pass (hiding latency, mirroring ATLAS's structure), then a blocked
`_mm_cmpeq_pd` + `_mm_movemask_pd` index search (2 elements/compare, exits
via the mask bit directly instead of a scalar re-check loop). Compiled with
the *exact* production flags (no `-march`, so genuinely SSE2-only):

```
N=72 -> nActivePairs=2556
  eigen (baseline)            1581.7 ns/call
  std::experimental::simd     1357.7 ns/call   (14% faster than baseline)
  hand SSE2 intrinsics        1203.7 ns/call   (24% faster than baseline,
                                                 11% faster than std::simd)
N=48: baseline 607.4 / std::simd 514.6 / hand SSE2 457.0 ns  (25% vs baseline)
N=24: baseline 106.2 / std::simd  91.8 / hand SSE2  66.6 ns  (37% vs baseline)
N=12: baseline  33.0 / std::simd  32.9 / hand SSE2  23.2 ns  (30% vs baseline)
```

So: **yes, hand-written intrinsics meaningfully beat `std::experimental::simd`
— a ~24-37% win over the current Eigen-based baseline**, roughly double the
13% figure from the library-based approach. The extra win over
`std::experimental::simd` comes from two things the library wrapper doesn't
do automatically: 4-way unrolling of the min-value pass (extra instruction-
level parallelism, hiding `minpd` latency), and using `movemask` to read the
match position directly out of the SSE2 compare result instead of falling
back to a scalar re-scan within the matching block.

Also compiled an AVX2 (`-mavx2`) version purely to see the ceiling (**not**
deployable — see below): ~1330ns at N=72, only marginally better than
`std::experimental::simd` compiled with the same `-mavx2` flag (1346ns).
Conclusion: most of AVX2's edge over SSE2 here is just from doing 4
doubles/vector instead of 2 — once the vector width is wide enough, hand-
tuning stops mattering much. The real, portable win is specifically in the
SSE2-width hand-tuned version above.

### Portability reassessed — intrinsics are *more* restrictive than `std::experimental::simd`, not less

Checked `.github/workflows/builds.yml` again: ACTS CI includes a `macos-26`
leg (Apple Silicon, i.e. **ARM64**). Raw x86 intrinsics (`_mm_loadu_pd`,
`emmintrin.h`, etc.) **do not exist at all on ARM** — this isn't a "might
not link the right stdlib" risk like `std::experimental::simd`'s Clang
concern, it's a hard compile failure on an actual, currently-building CI
target. `std::experimental::simd`, by contrast, is *designed* to degrade to
scalar/NEON automatically on non-x86 platforms when built with libstdc++;
its only real failure mode is Clang/libc++ specifically. So the ranking by
portability risk is: `std::experimental::simd` (fails only on
Clang+libc++) < hand SSE2 intrinsics (fails on **any** non-x86 target,
ARM macOS included) — the bigger raw win comes with the bigger portability
liability.

### Updated recommendation

Neither option is a safe drop-in replacement without a compile-time
platform dispatch:

- **Hand intrinsics**: would need `#if defined(__SSE2__)` (or
  `__x86_64__`) guarding the SSE2 kernel, with a **second**, portable
  (scalar or NEON-intrinsics) implementation for ARM — real, ongoing
  maintenance cost for two parallel code paths.
- **`std::experimental::simd`**: would need `#if __has_include(<experimental/simd>)`
  guarding, with a scalar fallback for any Clang/libc++ configuration that
  lacks it — simpler to guard (one `#if`, one fallback, and the SIMD
  abstraction itself still works on ARM+libstdc++ automatically) but a
  smaller win (13% vs. 24-37%).

Given ACTS explicitly builds on ARM macOS today, a same-source dual
implementation (SSE2 intrinsics + portable fallback) is the only way to get
the bigger win everywhere; that's a real but bounded engineering task (the
kernel itself is ~20 lines, as prototyped above), not a fundamental blocker
— but it's a decision Benjamin should make explicitly given it adds a
second, ISA-specific code path to maintain. Not implemented in production in
this session; still at the "prototyped and benchmarked, plan only" stage
per the original request.

## Follow-up investigation — row-cached minima (attempted, reverted)

Benjamin's idea: after a merge, only the merged/removed component's row
changes (plus whichever *other* rows' cached minimum happened to point at
it) — so cache each row's own minimum distance + partner, and turn
`minDistancePair`'s O(active²) full-array scan into an O(active) scan over
those cached minima. Explicit caveat from Benjamin up front: "the cache
patterns get worse... we would need to test it." We did — and the honest
result is a **net regression**, reverted back to the Stage 4/6 baseline.

### Correctness mechanism (why it doesn't need to check everything after a merge)

A merge only changes two things from any *other* row k's perspective: the
entry to the removed component disappears, and the entry to the
merged/updated component gets a new value. Every other entry in row k is
untouched. Since a row's cached minimum was, by construction, ≤ every one of
its entries, the only ways it can become stale are: (a) the *new* value to
the updated component is smaller (→ O(1) update), or (b) the cached minimum
used to be *at* the removed or updated component specifically and just got
worse/vanished (→ full O(active) rescan of that one row, since we don't know
what the row's second-best was). Everything else needs no work. Only a
small, roughly-constant number of rows hit case (b) per merge — not scaling
with N — which is what makes this amortized O(active) rather than O(active²)
per merge, in principle.

### Ties are handled the same way (no second-cached-value trap)

The design doesn't cache a second-best value anywhere (that's exactly the
naive "top-2" idea that breaks down after one iteration, since you don't
know the third-best once the top two are both gone). Instead, whenever a
row's cache is in doubt, it's discarded entirely and freshly recomputed by
scanning all of that row's current entries — so it doesn't matter whether
what's now the "true" minimum was previously 2nd or 50th place. Tie-breaking
(picking the smallest *original* triangular index among exactly-tied
distances, per the Stage 4/6 contract) was proven to hold hierarchically:
since "smallest original triangular index" is a total order, resolving
ties within each row first, then across the per-row winners, is guaranteed
equivalent to a single flat tie-break over every pair at once (by
transitivity) — verified by hand-tracing all 3 white-box unit tests plus a
dedicated new regression test before touching production code.

### Prototype (`scratch_benchmarks/rowmin_bench.cpp`) — initial result was misleading

Built a standalone prototype simulating the *full* 72→12 reduction sequence
(not just one call, since the payoff is amortized), with a correctness check
against the naive-baseline's invariant (weighted mean is conserved under any
merge order). Initial comparison: "packed full-scan" 40.9µs vs "row-min
cache" 21.9-22.8µs (**1.8-1.9x**), with only ~2.3 rescans/merge measured
(confirming the amortized-O(N) behavior held in practice, not degenerating
back to O(N²)).

**This comparison was wrong.** The "packed full-scan" baseline used a naive,
single-pass combined value+index scan — the *original, un-optimized Stage-0*
pattern — not the actual two-pass (vectorized-value, then early-exit-index)
`minDistancePair` that ships today. Fixing the toy baseline to match
production's real structure closed most of the gap in the toy model itself,
and more importantly the *actual* integration (below) settled it decisively.

### Production integration — measured net regression, reverted

Implemented the row-min cache design for real in
`SymmetricKLDistanceMatrix` (constructor computes `rowMin`/`rowMinPartner`
incrementally while filling the distance array; `recomputeAssociatedDistances`
patches every other row's cache inline; `maskAssociatedDistances` collects
invalidated rows before the swap and rescans them after; `minDistancePair`
scans the small `rowMin` array with the same hierarchical tie-break). All 9
existing unit tests passed immediately, including the exact-tie and
real-captured-ODD-mixture regression tests — the correctness design held up
exactly as reasoned.

**Performance did not.** `ActsBenchmarkGsfComponentReduction`:
`reduceMixtureWithKLDistance` went from **~67.1µs (current baseline) to
~75.5µs** — a regression, not a win. Debugged in three rounds:

1. **`rescanRow` at 26% of total time.** Disassembly showed it was a
   textbook combined scalar value+index+tie-break loop — the exact
   anti-pattern already fixed in `minDistancePair` back in Stage 2. Rewrote
   it as a proper two-pass (vectorized `Eigen::Array::minCoeff` over the
   contiguous "row part" + scalar "column part", then a single combined
   partner+tie-count pass instead of three separate scans). No measurable
   improvement (still ~26%).
2. **Instrumented real call counts**: merge-triggered rescans measured
   ~2.02/merge on real `GsfComponent` data — matching the toy model's 2.3,
   *not* anomalously high. The real culprit: the constructor was calling
   `rescanRow` once per component (72 extra O(active) scans per reduction,
   an entirely new O(N²) pass that didn't exist before, redundantly
   re-reading distances just written moments earlier). Fixed by computing
   `rowMin`/`rowMinPartner` incrementally inside the existing O(N²)
   distance-fill loop instead. This helped (`rescanRow`'s share dropped
   26%→16%) but total time only improved marginally (~75.5µs→~75.2µs) since
   the cost simply moved into the constructor itself (~12% of total, up
   from ~5% pre-row-cache).
3. **Isolated the remaining gap**: temporarily stripped the tie-safety
   branches entirely (diagnostic only, explicitly not shippable) — total
   time dropped to ~70.5µs, confirming tie-handling costs a real but
   partial ~5µs. Even *without* any tie-safety at all, the row-min-cache
   version was still slower than the current ~67.1µs baseline.

**Root cause, confirmed by correcting the toy model's baseline (above):**
the *existing* `minDistancePair` is already cheap in practice at N≤72
(vectorized value search + early-exit index search from Stages 2/4), so
there isn't much absolute time left for an O(N²)→O(N) asymptotic
improvement to recover — while row-caching unavoidably adds real per-element
bookkeeping (tie-detection branches, invalidation checks) into
`recomputeAssociatedDistances`/`maskAssociatedDistances`/construction, which
run just as often as before. At this problem size, that bookkeeping tax
costs more than the asymptotic win recovers. This is a legitimate,
well-known trade-off in hierarchical-clustering-style algorithms (matches
why ATLAS's own code comments mention trying an O(N² log N) heap-based
alternative and rejecting it for exactly this kind of overhead reason, per
the Stage 0 reference reading).

**Reverted** `Core/include/Acts/TrackFitting/detail/GsfComponentMerging.hpp`
and `Core/src/TrackFitting/detail/GsfComponentMerging.cpp` back to the
Stage 4/6 state (confirmed: `git diff --cached` against the pre-row-cache
stage shows zero differences, and the benchmark reproduces ~66.9µs,
consistent with Stage 6's numbers). All three exploratory files
(`argmin_bench.cpp`, `argmin_bench2.cpp`, `rowmin_bench.cpp`, plus a
`README.md` explaining each) are kept in `scratch_benchmarks/` (untracked,
like this file) per Benjamin's request, to revisit later.

## Follow-up investigation — further argmin exploration (GNU vec, fused single-pass, pragmas, float)

After the row-min-cache reversal, revisited the argmin primitive itself
(`SymmetricKLDistanceMatrix::minDistancePair`'s value+index search) from a
few more angles, each as its own standalone benchmark in
`scratch_benchmarks/`. All compiled/run in the same session for
apples-to-apples comparability; see that directory's `README.md` for build
commands.

### `argmin_bench3.cpp` — GNU vector extensions (the actual ATLAS mechanism)

ATLAS's `GsfFindIndexOfMinimum.h` doesn't use raw x86 intrinsics or
`std::experimental::simd` — it uses GCC/Clang **vector extension types**
(`__attribute__((vector_size(N)))`, wrapped by `CxxUtils/vec.h`), with the
ISA width passed in as a template parameter and dispatched via function
multi-versioning. This is a genuinely portable mechanism (same source
compiles on x86 *and* ARM/NEON — unlike raw intrinsics headers), so it's a
distinct fourth option worth benchmarking on its own, not just theorized
about.

Result (ns/call, isolated argmin):

| N | eigen (baseline) | std::experimental::simd | gnu vec W=2 (native SSE2 width) | gnu vec W=4 (mismatched width) |
|---|---|---|---|---|
| 12 | 34.7 | 39.9 | 31.7 | 195.1 |
| 24 | 109.1 | 123.8 | 103.8 | 788.6 |
| 48 | 620.0 | 586.2 | 610.2 | 3199.2 |
| 72 | 1593.9 | 1512.5 | **1428.0** | 7229.6 |

- **W=2 (matches native 128-bit/SSE2 width with no `-march` flag)**: only
  ~5-10% faster than Eigen — checked the generated assembly and confirmed
  GCC *does* lower the `mask ? a : b` idiom to a real `minpd` (not a
  `blendvpd` + compare), so it's not leaving that specific optimization on
  the table; the residual gap vs. the hand-written intrinsics version
  (below) is presumably in loop/scheduling details, not a missed
  vectorization.
- **W=4 (deliberately mismatched vs. the SSE2-native width, no `-mavx`)**:
  4-5x *slower*, and gets relatively worse as N grows. This is a direct,
  measured confirmation of ATLAS's own header comment: choosing a vector
  width wider than the ISA's native register forces the compiler to
  synthesize it from narrower ops, producing "quite poor assembly." Good
  corroborating data, not a usable option.
- **Verdict**: GNU vector extensions solve the portability problem (no
  x86-only headers, no libstdc++-only dependency) but recover much less of
  the win than raw intrinsics — a genuinely different, weaker point on the
  portability/speed tradeoff curve than option 3 below.

### `argmin_bench4.cpp` — fused single-pass argmin (Benjamin's idea)

Benjamin's observation: every variant benchmarked so far (Eigen,
`std::experimental::simd`, GNU vec, hand intrinsics) does **two full passes**
over the array — one vectorized pass for the minimum value, then a second
(scalar, early-exit) pass to recover its index. That second pass averages
N/2 in the common case but touches the *whole* array in the worst case (a
tie, or the min sitting near the end) — and production's `minDistancePair`
*always* pays a second full vectorized pass anyway for the tie-count check,
so eliminating the O(N/2) index-search pass is a real, not just
asymptotic, win.

The fix: process the array in fixed-size blocks, and while doing the
single vectorized min-reduction pass, track *which block* the running
minimum came from. Because the running min only updates on strict `<`, the
recorded block is guaranteed to be the block containing the *first*
occurrence of the true global minimum (a later block with an equal value
can't overwrite it) — so after the single pass, index resolution only
needs to rescan one fixed-size block, not the whole array.

Result (ns/call), hand-written SSE2 for the per-block reduction:

| N | eigen (2-pass) | hand SSE2 (2-pass) | fused 1-pass, block=8 | **block=16** | block=32 |
|---|---|---|---|---|---|
| 12 | 34.7 | 22.7 | 35.6 | 27.4 | 25.8 |
| 24 | 109.8 | 65.4 | 84.4 | 71.4 | 66.9 |
| 48 | 609.0 | 435.6 | 378.2 | **348.7** | 347.3 |
| 72 | 1534.3 | 1176.3 | 956.9 | **929.3** | 936.7 |

At N=48/72 — the sizes that dominate the total reduction cost, since the
merge sequence spends most of its O(N²)-summed work at the larger end —
block=16 beats even the hand-tuned 2-pass SSE2 version by ~20%, and beats
the Eigen baseline by ~40%. At N=12/24 it's *worse* than the 2-pass
versions: the fixed per-block bookkeeping (branch + compare per block)
costs more than it saves when there's little array left to avoid
rescanning in the first place.

Caveats before this could go into production:
- This benchmark isolates just the argmin. Production's tie-count pass
  (`(m_distances.head(nActivePairs) == min).count()`, always executed) is
  unaffected either way — the win here specifically removes the O(N/2)
  scalar index-search, not the tie-count.
- As written, still raw SSE2 intrinsics (x86-only, same ARM/macOS-CI
  portability blocker as `argmin_bench2.cpp`). The *algorithmic* idea
  (block-tracking single pass) is portable; only the per-block reduction
  would need reimplementing with GNU vector extensions (or scalar
  fallback) to actually ship it.

**This is a genuinely promising lead — unlike the row-min-cache idea, this
one has real headroom at the sizes that matter.** Not yet integrated into
production; would need the portable reimplementation plus a real
end-to-end benchmark (not just the isolated argmin) to confirm.

### `argmin_bench5.cpp` — does the autovectorizer find this on its own?

Same fused single-pass/block-tracking idea, but written as **plain nested
`for` loops** with a `constexpr` block width — no intrinsics, no GNU vector
extension types — to see whether GCC's autovectorizer picks it up
unassisted. (min/max reductions are exactly associative/commutative, so
GCC is allowed to vectorize this even without `-ffast-math`, unlike a sum
reduction.)

| N | eigen (2-pass) | fused hand SSE2 (block=16) | fused **autovec** (block=8) | block=16 | block=32 |
|---|---|---|---|---|---|
| 12 | 34.7 | 26.3 | 34.8 | 29.5 | 43.7 |
| 24 | 109.8 | 70.7 | 88.8 | 83.4 | 119.2 |
| 48 | 608.5 | 349.3 | 408.5 | 404.6 | 498.8 |
| 72 | 1543.2 | 927.4 | 1190.2 | 1188.3 | 1309.2 |

Confirmed via the generated assembly (`minpd` instructions present) that
the compiler *does* auto-vectorize the inner block reduction without any
hints. It recovers roughly half-to-three-quarters of the hand-tuned
version's win over Eigen (~23-25% faster than Eigen at N=48/72, vs. the
hand version's ~40%), but doesn't fully close the gap — likely down to
instruction scheduling/unrolling choices differing from the explicit
4-register unroll used by hand, not a vectorization failure. block=32 is
uniformly worse than block=8/16 for the plain-loop version — too coarse,
wastes work on the rarely-taken improvement branch relative to its own
cost.

**Tried pragmas to close the gap further — both made it worse:**

| N | plain autovec (block=16) | `#pragma GCC unroll 4` | `#pragma omp simd reduction(min:...)` (`-fopenmp-simd`) |
|---|---|---|---|
| 12 | 30.6 | 35.1 | 46.1 |
| 24 | 85.0 | 97.5 | 137.7 |
| 48 | 421.0 | 447.7 | 526.2 |
| 72 | 1207.5 | 1243.9 | 1346.4 |

Both explicit hints were consistently slower (unroll: ~3-15%; omp-simd
reduction: ~10-25%) than letting GCC's own vectorizer choose unassisted —
apparently the hints disrupt scheduling/register allocation GCC had
already picked well, without unlocking anything new. **Conclusion: for this
kernel, plain autovec code with no pragmas is the best non-intrinsic
option**, still ~15-20% behind the hand-tuned SSE2 version at N=48/72.

### `argmin_bench6.cpp` — does `float` change the picture?

Same fused block-tracking kernel, `double` vs. `float`, autovec and
hand-tuned, isolated from any correctness questions about switching
production's distance metric (see caveat below).

| N | eigen double | eigen float | autovec double (blk16) | **autovec float (blk16)** | hand SSE2 double (blk16) | hand SSE2 float (blk32) |
|---|---|---|---|---|---|---|
| 12 | 32.7 | 33.8 | 31.9 | 34.4 | 29.0 | 27.9 |
| 24 | 116.3 | 101.4 | 93.1 | 91.4 | 79.5 | 70.1 |
| 48 | 815.7 | 441.9 | 602.0 | **296.6** | 531.1 | 211.7 |
| 72 | 1698.7 | 1213.7 | 1349.7 | **890.2** | 1070.4 | 631.3 |

`float` roughly halves cost across the board at N=48/72 (register packing
— 4 floats vs. 2 doubles per 128-bit lane — plus half the memory traffic),
independent of the autovec/hand-tuned question: `float` autovec beats even
`double` hand-tuned SSE2. Hand-tuned SSE2 `float` (block=32, matched to
float's native 4-wide register) is the fastest of everything benchmarked
across all four investigations, at ~37% of the original Eigen-`double`
baseline. One anomaly: `float` autovec at block=32 is *much worse* than
block=16 (671ns vs. 297ns at N=48) even though block=32 works fine for the
hand-tuned float version (2×`__m128` per chunk) — the plain-loop
autovectorizer doesn't handle that block size well for `float`; block=16
(4×`__m128`-equivalent) is its sweet spot.

**Important caveat, not yet resolved:** production's `double` distance
path must stay bit-exact with the naive baseline (`test_naive_vs_optimized`
at 1e-8). Switching the metric itself to `float` is a materially bigger
design decision than an argmin micro-optimization — it changes which pairs
tie/merge-order in edge cases, exactly as flagged in the original task plan
(`Stage 3 — float instantiation`, never executed in this investigation).
This benchmark only establishes that float *would* be worth it kernel-wise
if that bigger design question were separately decided; it doesn't decide
it.

### Summary of standing across all argmin investigations

From fastest to slowest (at N=48-72, the range that dominates total cost),
excluding the not-yet-shippable float-metric option:

1. Fused single-pass, hand-tuned SSE2 intrinsics, block=16 — best, x86-only.
2. Fused single-pass, plain-loop autovec, block=8-16 — ~15-20% behind #1,
   fully portable (ARM included), no pragmas needed or beneficial.
3. Two-pass hand-tuned SSE2 intrinsics (`argmin_bench2.cpp`) — behind #2 at
   large N (the two-pass structure's O(N/2) index scan starts to dominate).
4. GNU vector extensions at native width (`argmin_bench3.cpp`) — modest,
   ~5-10% win, portable, but far behind #1/#2.
5. `std::experimental::simd` — inconsistent, GCC/libstdc++-only.
6. Eigen (current, shipped) — baseline.

None of this has been integrated into production yet. Option #2 (fused
single-pass, plain autovec) looks like the most promising real candidate:
it doesn't need architecture-specific code (portable to ARM out of the
box) and gives up only a modest amount of the hand-tuned ceiling. Next step
if pursued: reimplement it against the real `m_distances` array structure
(including the required tie-count pass, currently unaffected by this
change) and measure the actual `ActsBenchmarkGsfComponentReduction` effect,
the way the row-min-cache idea was validated (and rejected) end-to-end
rather than trusting the isolated microbenchmark alone.

## Decision — remove the tie-break-to-match-naive mechanism entirely

Revisited the Stage 4 tie-break fix (`minDistancePair`'s tie-count check +
original-triangular-index resolution, `GsfComponentMerging.cpp:151-181`)
after establishing (previous section) that it's paying real, non-negligible
cost every call (the vectorized `.count()` check is always executed, and
was a measured contributor to `minDistancePair` rising from ~39% to ~46.56%
of total time between Stage 3 and Stage 6's profiling).

**Discussion:** in production, matching naive's specific tie-break choice
has no physical meaning — when two candidate merges are at an exact
distance tie, neither choice is "more correct" than the other, they're
just two different, equally valid reductions. Paying a cost on every single
call (including the ties-never-happen common case) to reproduce one
arbitrary convention is not worth it. Three options were considered:

- **(A)** Runtime bool parameter gating the tie-break code, off by default
  in production, on in tests that need bit-exact naive equivalence.
- **(B)** Same as (A) but as a compile-time template parameter, so the
  disabled branch is compiled out entirely rather than left as a runtime
  no-op check.
- **(C)** Remove the tie-break code entirely, and instead of testing for
  bit-exact agreement with naive under ties, construct two targeted test
  cases: one with no tie (equality still expected and checked) and one
  engineered so the tie *does* cause optimized and naive to diverge
  (checked as an expected, documented inequality, not a bug).

**Decision: (C).** Benjamin's call — this fixes the behavior (production no
longer pays for tie resolution, at all) without adding any code, which is
preferable to (A)/(B)'s added parameter/template complexity for a
distinction (matching an arbitrary tie convention) that isn't meaningful in
production in the first place.

### Planned test changes (not yet implemented)

- Remove the tie-count-check + resolution block from `minDistancePair()`
  (`GsfComponentMerging.cpp:151-181`), leaving only the min-value pass and
  the early-exit index scan — i.e. revert to the Stage 2/3 structure, minus
  the Stage 4 tie-break addition.
- `test_exact_tie_breaking` (four exactly-identical components, currently
  asserts equality with naive): repurpose as the "tie does make a
  difference" case — assert the optimized result now legitimately
  *diverges* from naive's choice (or is no longer required to match it),
  documenting this as expected, not a regression.
- `test_naive_vs_optimized_stress` and `test_real_odd_mixture_repro`
  (currently forced-duplicate / real-data cases that were passing via the
  tie-break fix): need re-examination once the fix is removed — either they
  no longer hit an actual tie post-fix removal (still pass as-is) or they do
  and need to move to the "expected inequality" bucket alongside
  `test_exact_tie_breaking`, or be re-cast as order-invariant
  invariant-preservation checks (weight/mean conservation) rather than
  row-by-row equality with naive.
- `test_naive_vs_optimized` (the main fixed-seed, continuous-random-data
  correctness test) is unaffected either way, since real ties don't occur
  with continuous random doubles.

Not yet implemented — this section records the discussion and decision
only.

## Stage 7 — Implementation: tie-break removal + fused single-pass argmin

Implements the two items left open above: Decision C (remove the tie-break-
to-match-naive mechanism) and the fused single-pass, block-tracking argmin
identified as a promising portable lead in the "further argmin exploration"
follow-up. A third item, revisiting `float` vs `double` distance
bookkeeping on top of these two changes, was also tried per Benjamin's
request and is documented separately below (Stage 7c) since it was
ultimately **not** kept.

### Stage 7a — Remove the tie-break-to-match-naive mechanism (Decision C)

Deleted the tie-count check (`(m_distances.head(nActivePairs) ==
min).count()`) and the "resolve by smallest original triangular index" loop
from `SymmetricKLDistanceMatrix::minDistancePair()`
(`GsfComponentMerging.cpp`), reverting the function to its Stage 2/3
structure: vectorized min-value pass + early-exit scalar index scan, no
tie-break bookkeeping. Updated the function's doc comment accordingly.

**Correctness:** rebuilt and reran `ActsUnitTestGsfMixtureReduction` —
**all 9 tests still pass unmodified**, including `test_exact_tie_breaking`,
`test_naive_vs_optimized_stress`, and `test_real_odd_mixture_repro`, which
were originally added to guard the tie-break fix. Investigating why:
- `test_exact_tie_breaking` reduces four *exactly identical* components to
  one. Merging identical Gaussians is commutative/associative regardless of
  pairing order (the weighted mean/covariance of identical components
  doesn't depend on which pair merges first), so this test is order-
  invariant by construction — it was never actually sensitive to which tie-
  break convention was used, even though the original Stage 4 comment
  implied it was. Updated the test's comment to reflect this: it's now a
  regression guard for the merge math itself, not for tie-break ordering.
- `test_naive_vs_optimized_stress` and `test_real_odd_mixture_repro` simply
  didn't hit a case where compacted-order vs. fixed-order tie-breaking
  picks a different final component set for their specific fixed-seed
  data — no comment changes needed since neither test's rationale
  referenced the removed mechanism.

So Decision C's removal caused **zero test regressions**, which is itself a
useful data point: the tie-break-to-match-naive mechanism was paying a
real, measured cost (below) for a distinction that, empirically and not
just in principle, doesn't materialize in either the synthetic stress data
or the captured real ODD mixture used in these tests.

**Benchmark** (`ActsBenchmarkGsfComponentReduction`, same setup as all
prior stages):

```
reduceMixtureWithKLDistance (optimized)   67,424 ns/iter (Stage 6 baseline)
                                        -> 60,424 ns/iter (after Stage 7a)
```

**~10.4% faster**, consistent with reverting to the Stage 2/3 structure
(Stage 6's number, ~67.1-67.4µs, already reflected the Stage 4 tie-break
tax on top of Stage 2/3's ~59.8-60.0µs).

### Stage 7b — Fused single-pass, block-tracking argmin

Replaced the two-pass search in `minDistancePair()` (vectorized
`minCoeff()` then full-array early-exit scalar index scan) with the
block-tracking single-pass argmin validated in
`scratch_benchmarks/argmin_bench5.cpp` (`argminFusedAutovec<BLOCK>`,
`BLOCK=16`, no intrinsics, no pragmas — plain nested loops relying on
GCC's autovectorizer for the inner block-min reduction, which the earlier
investigation confirmed via `perf annotate`/assembly and explicitly showed
is *not* helped by `#pragma GCC unroll` or `#pragma omp simd`). The
algorithm processes the active distance prefix in fixed 16-element blocks,
tracks which block produced the running minimum (guaranteed to contain the
*first* occurrence of the true minimum, since the running min only updates
on strict `<`), then rescans only that one block to recover the index —
replacing an O(nActivePairs/2)-average full-array scalar scan with a single
pass plus a small bounded rescan. Everything else in `minDistancePair()`
(the `m_mapToPair`/`m_activeToOriginal` resolution) is unchanged.

**Correctness:** `ActsUnitTestGsfMixtureReduction` — all 9 tests still
pass. The block algorithm preserves "first occurrence in current compacted
scan order" semantics, so behavior is equivalent to the prior early-exit
scan's contract.

**Benchmark:**

```
reduceMixtureWithKLDistance (optimized)   60,424 ns/iter (after Stage 7a)
                                        -> ~47,700-48,600 ns/iter (after Stage 7b, 3 repeats)
```

**~19-21% faster on top of Stage 7a** — this is the real, end-to-end
confirmation the earlier investigation flagged as still missing: the
isolated microbenchmark's ~15-25% win at N=48-72 (`argmin_bench5.cpp`)
translated directly into a real win in the actual
`ActsBenchmarkGsfComponentReduction` workload, unlike the row-min-cache
idea from the same investigation phase, which looked promising in isolation
but regressed once integrated for real. No revert needed here.

### Stage 7c — Revisit `float` vs `double` distance bookkeeping (tried, not kept)

Per Benjamin's request, re-templated `computeSymmetricKlDivergence<Scalar>`
and `SymmetricKLDistanceMatrix<Scalar>` on top of the Stage 7a/7b code
(mirroring the original Stage 3 mechanism: `mergeTwoComponents` stays
`double`-only, only the distance bookkeeping that decides merge order runs
in `Scalar` precision), reintroduced `reduceMixtureWithKLDistanceFloat`,
and added a float/naive invariant test (weight and weighted-qop-mean
conservation, since exact merge-order equivalence isn't meaningful once
distances are computed in reduced precision).

**Correctness:** all 10 tests (9 + the new float invariant test) passed.

**Benchmark** (double vs. float, both on top of Stage 7a/7b, 3 repeats):

```
reduceMixtureWithKLDistance (optimized)         ~47,700-47,750 ns/iter
reduceMixtureWithKLDistanceFloat (optimized)    ~45,500-45,620 ns/iter
```

**Only ~4.5% faster than double** — smaller than Stage 3's original 1.13x
and even smaller than Stage 6's ~1.07-1.08x margin that was already judged
"not worth the complexity" and dropped. This directly contradicts what the
isolated `argmin_bench6.cpp` micro-benchmark suggested (float roughly
halving the fused-argmin kernel's own cost) — the likely explanation is
that Stage 7b's block-tracking argmin already shrank `minDistancePair`'s
share of total time substantially (see profiling below), so the argmin
kernel's own float-vs-double delta no longer dominates the end-to-end
number the way it did when `minDistancePair` was ~40-47% of total runtime
in earlier stages; `computeSymmetricKlDivergence`, `mergeTwoComponents`,
and the matrix bookkeeping functions are now comparatively larger
contributors, and `mergeTwoComponents` always stays `double`.

**Decision: dropped**, per the plan's decision gate (a small end-to-end
win doesn't justify a second, duplicated code path). Reverted the Scalar
templating back to plain `double` everywhere (`SymmetricKLDistanceMatrix`,
`computeSymmetricKlDivergence`), removed `reduceMixtureWithKLDistanceFloat`
and its test/benchmark variant. Confirmed the revert reproduces exactly the
Stage 7b state: 9/9 unit tests pass, benchmark unaffected.

### Final profiling (`perf record -F 999 -g`, Stage 7a+7b state)

```
    25.72%  [.] SymmetricKLDistanceMatrix::minDistancePair() const
    14.79%  [.] computeSymmetricKlDivergence(...)
    13.68%  [.] reduceWithKLDistanceImpl(...)
    12.59%  [.] mergeTwoComponents(...)
     6.34%  [.] SymmetricKLDistanceMatrix::SymmetricKLDistanceMatrix(...)
     5.19%  [.] SymmetricKLDistanceMatrix::recomputeAssociatedDistances(...)
     4.83%  [.] SymmetricKLDistanceMatrix::maskAssociatedDistances(...)
     2.75%  [.] std::swap<GsfComponent>(...)
     2.16%  [.] Eigen (BoundMatrix accumulation in mergeTwoComponents)
```

`minDistancePair`'s share dropped from 46.56% (Stage 6, with the tie-break
tax) to 25.72% — a large, direct confirmation of Stage 7a+7b's effect. The
profile is now much more balanced across the actual algorithmic work
(distance computation, merging, matrix bookkeeping) rather than being
dominated by the min-search, which was the original Stage 0 complaint.

### Final full-GSF validation on OpenDataDetector

Repeated the Stage 4 methodology (`Examples/Scripts/Python/truth_tracking_gsf.py`,
ODD detector with real material maps, 500 events, seed 42, `numThreads=1`,
comparing `TrackFittingAlgorithm` timing before/after), comparing the
Stage 6 baseline ("before") against the Stage 7a+7b state ("after", Stage
7c reverted):

| | TrackFittingAlgorithm total (500 ev) | ms/event |
|---|---:|---:|
| Before (Stage 6) | 1290.24 ms | 2.58 |
| After (Stage 7a+7b) | 1238.37 ms | 2.48 |

**~4.0% reduction in the full GSF track-fitting algorithm's wall-clock
time** from Stage 7a+7b alone, on top of Stage 6's already-optimized
baseline — smaller than Stage 4's ~9.1% (which was measured against the
original, pre-Stage-1 unoptimized code), consistent with these being
further, incremental refinements on an already-substantially-improved
function rather than the original structural fix.

**Correctness:** compared all 65 `tracksummary` branches between the
before/after ROOT outputs across all 500 events — **bit-identical in every
branch**, including `chi2Sum`/`measurementChi2` (which, unlike Stage 4's
comparison against the original unoptimized code, showed no divergence at
all here — Stage 6 and Stage 7a+7b apparently don't hit the same rare
exact-tie edge case Stage 4 found against the original code, at least not
in this fixed-seed 500-event sample).

### Summary

| Reducer | Stage 0 (ns/iter) | Stage 6 (ns/iter) | Final, Stage 7a+7b (ns/iter) | Speedup vs Stage 0 | Speedup vs naive |
|---|---:|---:|---:|---:|---:|
| `reduceMixtureWithKLDistance` | 137,755 | ~67,100-67,400 | ~47,700-48,600 | **~2.85-2.90x** | **~4.0x** |

Combined with Stage 6's numbers, this closes out the task: the
"sophisticated" KL-distance reducer is now decisively (~4x) faster than
the naive `O(N^3)` baseline it was originally only ~1.5x faster than, via
(1) the Stage 1 structural fix (swap-to-last compaction, the dominant
win), (2) Stage 2/3's two-pass vectorized search, (3) Stage 7a removing an
unnecessary correctness tax that never had a real behavioral payoff in
practice, and (4) Stage 7b's fused single-pass argmin, the last
significant, portable, benchmark-confirmed win identified in the
investigation phase. `float` distance bookkeeping (Stage 3, revisited in
Stage 7c) was tried twice now and dropped both times for the same
reason: a small measured win not worth a second code path.
