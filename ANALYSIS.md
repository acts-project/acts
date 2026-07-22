# Gen3 navigation performance — profiling and fixes

Gen3 geometry propagation through the Open Data Detector was ~14 % slower than the equivalent
Gen1 run, almost entirely in navigation. This note records the profiling that localized the gap
and the changes that close most of it. Net effect (100k single-threaded ODD propagations): gen3
went from **99.3 → 15.2 heap allocations per propagation** and from **~4.55 s → ~4.10 s**, versus
gen1 at 13.2 alloc/prop and ~4.02 s — i.e. from ~14 % slower to ~2 %, with Gen1 unchanged.

Changes, each a separate commit and independently measured:

1. **Remove per-propagation heap allocations from Gen3 navigation** — allocation-free candidate
   de-duplication in `NavigationStream::initialize`, and small-buffer/inline storage for the
   type-erased navigation policy states (`std::any` → `Acts::AnyBase`). 99.3 → 15.2 alloc/prop.
2. **Store `MultiNavigationPolicy` child states linearly** — drop the vector-of-handles in favour
   of index addressing; removes a hack and shrinks the state.
3. **Skip the per-step `isValid()` when the policy state is the default sentinel** — the always-
   true validity check was pure per-step overhead for policies that don't override it. −2.4 %
   instructions.

A byte-arena rewrite of the state manager was prototyped and **reverted** (measured slightly
worse — see the follow-up section). All changes are Gen3-only; Gen1 is byte-for-byte unchanged and
the navigation unit tests pass.

## Setup

- Worktree: `.claude/worktrees/gen3-prop-profiling` (branch off `main` @ `081ddad192`).
- Lean **RelWithDebInfo** build, gcc-15 (matches spack deps ABI), python 3.14, in-worktree
  `build/` (Python bindings + DD4hep + ODD; assertions OFF for realistic perf).
- Env: `source wenv.sh` (gcc-toolset-15 + spack `ci-dependencies` + build + venv, `ODD_PATH` set).
- Driver: `prof_propagation.py --mode {gen1,gen3}` — ODD detector, muon particle gun
  (η∈[-4,4], pT 1–100 GeV, 1000/ev), `EigenStepper`, constant 2 T field, single thread,
  **IO disabled** (sterile logger, no writers), **material off** (isolates geometry/navigation;
  also avoids the ROOT-plugin dependency). Same stepper/field/gun for both modes — the only
  difference is `getOpenDataDetector(gen3=…)`, i.e. the geometry + navigator dispatch.

## Headline result

100k propagations (100 ev × 1000 trk), mean of 3 runs:

| mode | propagation run | geometry build |
|------|-----------------|----------------|
| gen1 | **3.98 s**      | 0.62 s |
| gen3 | **4.55 s**      | 0.61 s |

**Gen3 is ~14–15 % slower** in propagation (geometry build time is equal). Reproduced at
200k (8.13 s vs 9.35 s) with perf profiles below.

## Flat profile (self-time, `perf -F 4000`, 200k propagations)

Shared cost (identical in both, ~40 %): `EigenStepper::step`, `CylinderSurface::intersect` +
`intersectionSolver`, `PlanarHelper::intersect`, `PlaneSurface/DiscSurface::intersect`,
`Surface::localToGlobalTransform`, `Eigen::Transform::inverse`, `updateSingleSurfaceStatus`,
`checkPathLength`. These are stepping + surface intersection and are not the gen1/gen3 delta.

The delta is entirely in **navigation** and **heap churn**:

| Gen1 navigation | % | Gen3 navigation | % |
|---|---|---|---|
| `Layer::compatibleSurfaces` (+lambda) | 9.3 | `NavigationStream::initialize` | 5.5 |
| `TrackingVolume::compatibleLayers` | 1.4 | `CylinderNavigationPolicy::initializeCandidates` | 1.6 |
| `getNextTargetGen1` | 1.1 | `getNextTargetGen3` | 1.3 |
| `resolveSurfaces` | 0.8 | `handleSurfaceReached` | 1.2 |
| SurfaceGridLookup | 0.9 | `SurfaceArrayNavigationPolicy::initializeCandidates` | 1.1 |
| | | `resolveCandidates` | 1.0 |
| | | candidate `insertion_sort` | 0.9 |
| | | `MultiNavigationPolicy::createState` | 0.8 |
| | | SurfaceGridLookup | 1.2 |
| **≈ 13.5 %** | | **≈ 16–17 %** | |
| heap `malloc/free/memmove` ≈ **0.7 %** | | heap `malloc/free/memmove` ≈ **4.5 %** |

So gen3 spends ~3 % more in navigation logic **and ~4 % more in the allocator** — together
consistent with the ~14 % wall-clock gap.

## Root cause (verified against source)

`Navigator::getNextTargetGen3` (`Core/src/Navigation/Navigator.cpp:501`) is efficient in the
common case: while the current volume's policy stays valid it just advances an index
(`++navCandidateIndex`) and only calls the expensive `resolveCandidates` on volume
entry / candidate exhaustion — **not every step**. The cost is in what `resolveCandidates`
(`Navigator.cpp:562`) → `NavigationStream::initialize` (`Core/src/Navigation/NavigationStream.cpp:20`)
does each time it *is* called, and the allocations around it:

1. **`NavigationStream::initialize` does two sorts + two dedup passes per rebuild**
   (`NavigationStream.cpp:30` stable_sort by surface ptr, `:34` unique/erase, `:96` sort by
   path length, `:102` unique/erase again, `:115` resize). Gen1's `Layer::compatibleSurfaces`
   sorts once into a stack `small_vector<NavigationTarget,10>`.

2. **Heap allocations on the rebuild path** (the ~4.5 % allocator time):
   - `std::vector<NavigationTarget> additionalCandidates` allocated fresh inside
     `initialize` (`NavigationStream.cpp:42`) — only populated for surfaces with two valid
     intersections, but constructed unconditionally.
   - `SurfaceArrayNavigationPolicy::initializeCandidates` calls `SurfaceArray::neighbors()`,
     which returns a **heap `std::vector<const Surface*>`** every rebuild.
   - `MultiNavigationPolicy::createState` allocates a `std::vector<NavigationPolicyState>`
     and type-erased state per volume entry (`MultiNavigationPolicy.cpp:78`).
   - Copy-out: surviving candidates are copied `stream.candidates()` → `state.navCandidates`
     (`Navigator.cpp:645-651`); `navCandidates` is `small_vector<…,10>` and spills to heap
     when a volume yields >10 candidates (ODD sensitive layers routinely do).
   - `state.stream.candidates()` itself is `reserve(50)` once (`Navigator.cpp:107`) and
     `reset()` keeps capacity — good — so the stream vector is *not* the main churn source;
     the per-rebuild temporaries above are.

Gen1 avoids most of this: `compatibleSurfaces`/`compatibleLayers` return stack small_vectors
by value (RVO), one sort, no per-volume policy-state allocation, no copy-out stage.

### Empirical allocation count (LD_PRELOAD malloc counter)

Counting real `malloc` calls with `mallocount.so`, at two workloads per mode; the difference
isolates the per-propagation rate (subtracting fixed geometry-build + interpreter startup):

| mode | mallocs @100 prop | mallocs @100k prop | **per propagation** |
|------|-------------------|--------------------|---------------------|
| gen1 | 4,344,613         | 5,660,556          | **≈ 13.2** |
| gen3 | 4,876,526         | 14,800,760         | **≈ 99.3** |

**Gen3 issues ~7.5× more heap allocations per propagation (~86 extra each).** `free` counts
track `malloc` (14.80 M vs 14.75 M at 100k) so it is balanced churn, not a leak; `calloc`
(~2.5 k) and `realloc` (~21 k) are startup-only and identical across modes. Individually these
allocations are cheap (mostly glibc tcache hits, hence "only" ~4.5 % allocator self-time in
the perf flat profile) — but 86 extra small allocations per propagation × the candidate-rebuild
temporaries above is exactly the churn that opportunity #1 removes. This is the single clearest
lever: eliminating the per-rebuild temporaries should recover most of the gen3-vs-gen1 gap.

## Changes implemented & measured

Two responsibly-local, gen3-only changes (no API used outside gen3 navigation was
touched; gen1 is provably unaffected). Verified with the navigation unit tests
(`ActsUnitTestNavigationStream` passes; `ActsUnitTestNavigationPolicy` passes except one
pre-existing case that needs `ACTS_ENABLE_LOG_FAILURE_THRESHOLD=ON`, unrelated to these edits).

1. **`NavigationStream::initialize` — drop the allocating `stable_sort`, reuse a scratch buffer**
   (`Core/src/Navigation/NavigationStream.cpp`, `NavigationStream.hpp`).
   The pre-intersection surface-pointer dedup used `std::ranges::stable_sort`, which allocates a
   `std::_Temporary_buffer` on every call (confirmed as a top allocator by backtrace
   attribution). Replaced with an in-place first-occurrence dedup scan that reproduces the exact
   previous "first-wins" result (important when a surface is added both by a policy and as an
   external/free surface with a different boundary tolerance) — no sort, no allocation. Also made
   the `additionalCandidates` temporary a reused member buffer.

2. **`MultiNavigationPolicy` — inline the per-volume policy-state list**
   (`Core/include/Acts/Navigation/MultiNavigationPolicy.hpp`, `.cpp`).
   `State::policyStates` and the `createState` temporary were `std::vector<NavigationPolicyState>`
   (heap for the typical 2-policy volume). Switched to
   `boost::container::small_vector<NavigationPolicyState, 4>` (the idiom already used in
   `Navigator.hpp`), removing two allocations per volume entry.

Measured (100 ev × 1000 trk = 100k propagations; malloc/prop via two-workload difference):

| | malloc / prop | run mean |
|------|--------------|----------|
| gen1 before / after | 13.17 / **13.17** (unchanged) | 4.017 s / **4.023 s** |
| gen3 before / after | 99.35 / **44.57**  (**−55 %**) | 4.547 s / **4.290 s** (**−5.6 %**) |

Heap allocations per gen3 propagation more than halved; wall-clock improved 5.6 %, shrinking the
gen3-vs-gen1 gap from ~13 % to ~6.6 %. gen1 is unchanged, confirming the edits are gen3-local.

3. **`NavigationPolicyStateManager` — small-buffer type erasure (hidden in the manager)**
   (`Core/include/Acts/Navigation/INavigationPolicy.hpp`).
   The manager stored states as `std::vector<std::any>`; every `std::any::emplace` of a state
   larger than the ~16-byte `std::any` SBO heap-allocated (the 88-byte
   `MultiNavigationPolicy::State` did so on every volume entry — the dominant remaining allocator
   per backtrace attribution). Switched the storage to `std::vector<Acts::AnyBase<128>>` (ACTS's
   configurable-SBO any); 128 B holds every policy state inline (largest is 88 B), so pushes no
   longer allocate. A `static_assert(sizeof(T) <= kNavigationPolicyStateSbo)` in `pushState`
   makes a future oversized state a compile error rather than a silent heap fallback. The stack
   vector reserves its typical depth lazily on first push, so the Gen1 navigator — which carries
   the manager but never pushes — still allocates nothing. **No navigation policy changed**: the
   type erasure stays hidden behind `NavigationPolicyState::as<T>()` / `pushState<T>()`.

Final measured state (100k propagations):

| | malloc / prop | run mean |
|------|--------------|----------|
| gen1 baseline / final | 13.17 / **13.17** (unchanged) | 4.017 s / 4.033 s |
| gen3 baseline | 99.35 | 4.547 s |
| gen3 after #1+#2 (dedup + small_vector) | 44.57 | 4.290 s |
| gen3 after #3 (SBO state manager) | **15.17** | **≈ 4.20 s** |

Net: gen3 heap allocations **−85 %** (99.4 → 15.2 / prop, now within ~2 of gen1), wall-clock
**≈ −7 %**, shrinking the gen3-vs-gen1 gap from ~13 % to ~4 %. gen1 is byte-for-byte unchanged.
Validated by `ActsUnitTestNavigationStream`, `ActsUnitTestNavigationPolicyState`, and
`ActsUnitTestNavigationPolicy` (all pass; the one skipped case needs
`ACTS_ENABLE_LOG_FAILURE_THRESHOLD=ON`, unrelated).

## Follow-up: state-manager call overhead (ceiling, linear-push, arena)

After the allocation work, gen3 sits at 15.2 malloc/prop / ~4.2 s vs gen1 13.2 / ~4.0 s.

**Ceiling.** Temporarily bypassing the entire policy-state subsystem (nothing reads state
contents yet) gets gen3 to **3.99 s / 14.2 malloc/prop** — parity with gen1. So the whole
remaining prize is ~5 % CPU + ~1 malloc/prop, and it is *per-call* overhead, not heap traffic.
Profiling attributes it to `MultiNavigationPolicy::createState` (2.15 %, per volume entry) and
`isValid` (1.37 %, per step) plus `initializeCandidates`/`popState`.

**Linear-push (committed).** `MultiNavigationPolicy::State` no longer stores a vector of child
handles; children are pushed contiguously and addressed by index via new
`NavigationPolicyState::index()`/`atIndex()`. Removes the hack, shrinks the state 88→4 B, makes
it trivially destructible, and drops ~0.4 % of gen3 instructions (40.24e9 → 40.08e9). Dispatch is
bit-identical; all navigation unit tests pass. (The child *count* stored here was later dropped
too — see the sentinel change — since it always equals the number of contained policies.)

**Full byte-arena (implemented, measured, reverted — negative result).** Replacing
`vector<AnyBase<128>>` with an offset-based byte arena (placement-new into one reused buffer +
`{offset,dtor,typeId}` frame table, type-checked reinterpret access) was **~0.3 % *more*
instructions** (40.20e9 vs 40.08e9) with identical allocation. The ceiling analysis explains why:
the cost is per-call overhead (`createState`/`isValid`), not the storage mechanism, so a
storage refactor cannot capture it — and the arena trades `AnyBase`'s cheap one-time handler
indirection for `resize`/bounds-check/alignment work plus real fragility (raw-byte lifetimes,
trivial-relocation assumption, no-pointer-across-push). Reverted; kept the simpler `AnyBase`
version.

**Sentinel isValid-skip (committed).** The remaining per-step cost was `isValid`, which the
navigator calls on *every* step to check the resolved candidates are still valid — but no leaf
policy overrides it (the default just returns `true`), so for the ODD it is always-true dead
weight. Rather than speed up the dispatch, skip it: a policy with a meaningful `isValid` must push
a real state to check against, whereas trivial policies push the default `EmptyState`. So the
*state itself* is the signal — `NavigationPolicyState::isDefault()` (true for `EmptyState`/no
state) tells the navigator when it can skip `isValid` entirely. This is footgun-free (a real
`isValid` ⟹ a real state ⟹ auto-detected; nothing to remember to set), verified by
`ConeValidityPolicy` (which implements `isValid`) still forcing re-resolutions unchanged.
`MultiNavigationPolicy` pushes `EmptyState` when all children are default and its marker `State`
only when one can invalidate; child states are addressed by policy count, so the stored
`childCount` is dropped. Measured: gen3 **40.12e9 → 39.17e9 instructions (−2.4 %)**, ~4.20 s →
~4.10 s; gen1 unchanged.

## Opportunities identified (status)

Items 1–2 and the per-step `isValid` (item 4's second half) are **done in this PR**; item 3 and
the shared-inverse work below are **follow-ups not included here**. Note: `SurfaceArray::neighbors`
turned out to already return a `std::span` into precomputed storage (no allocation), so only the
other two sub-points of item 1 were real.

1. **Kill the per-rebuild temporaries** (highest value, ~4 % allocator). ✓ done
   - `SurfaceArray::neighbors()` returning a fresh `std::vector` each call → add an overload
     that appends into the caller's `AppendOnlyNavigationStream` (or a reused scratch buffer
     in `Navigator::State`). This is the single biggest allocation on the path.
   - `additionalCandidates` in `NavigationStream::initialize` → make it a reused member/scratch
     buffer (`clear()` instead of construct), or avoid it by handling the second intersection
     in place; it is empty for the vast majority of surfaces.
   - `MultiNavigationPolicy` policy-state vector → reuse a buffer in the state manager instead
     of allocating per volume entry.

2. **Reduce the sort/dedup work in `NavigationStream::initialize`.**
   Two sorts + two `unique` passes over the candidate set each rebuild. The initial
   stable_sort-by-pointer + unique is a de-dup that could be folded into the single
   path-length sort (dedup after the sort you already need), removing one O(n log n) pass and
   one O(n) pass. `insertion_sort` shows up at 0.9 % — candidate counts are small, so a
   branch-lean single sort matters.

3. **Grow `state.navCandidates` inline capacity or navigate straight off the stream.**
   The `stream.candidates()` → `navCandidates` copy (`Navigator.cpp:645`) duplicates the
   candidate list every rebuild and spills to heap for busy layers. Either bump the
   `small_vector` inline size to cover typical ODD layers, or iterate the stream directly and
   drop the second container.

4. **`CylinderNavigationPolicy::initializeCandidates` (1.6 %)** is allocation-free already;
   its cost is geometric (portal intersections). Lower priority, but the `isValid` check runs
   every step — worth confirming it isn't recomputing intersections that the rebuild already has.

## Shared opportunities (help both gen1 and gen3, ~40 % of runtime)

- `CylinderSurface::intersect` + `intersectionSolver` (~9–11 %): pure geometry; already
  inverse-free.
- `Eigen::Transform<double,3,Affine>::inverse` (~1.8–2.6 %): this is **not** on the intersection
  path — the cylinder/plane/disc `intersect` methods already avoid the inverse (they do the
  boundary check with the rotation transpose, `tMatrix.block<3,2>(0,0).transpose()`, "built-in
  local to global for speed reasons"). The inverse comes from `Surface::globalToLocal`
  (`localToGlobalTransform(gctx).inverse() * position`), called once per **reached** surface from
  bound-state conversion (`TransformationHelpers.cpp:49`), the free→bound Jacobian
  (`JacobianEngine.cpp:133`, `Surface.cpp:49`), and material interaction
  (`PointwiseMaterialInteraction.cpp:42`, off in these runs). Call sites: `PlaneSurface.cpp:80`,
  `DiscSurface.cpp:95/131/207/243`, `CylinderSurface.cpp:138/162/340`.
  Two independent, orthogonal fixes (both help gen1 and gen3 equally — this is *not* part of the
  gen3 gap):
  1. **Cheaper inverse, no caching:** surface placements are rigid isometries, but `.inverse()`
     defaults to the general affine inverse. Call `.inverse(Eigen::Isometry)` (→ `Rᵀ`, `-Rᵀt`) —
     exact, much cheaper, purely local edits at the call sites.
  2. **Cache the inverse:** for statically-placed surfaces (`m_placement == nullptr`) the
     transform is constant, so store its inverse next to `m_transform` (computed at construction)
     and expose a `globalToLocalTransform(gctx)` accessor. Caveat: for **alignable** surfaces
     (`m_placement != nullptr`, transform depends on `GeometryContext`) the inverse must be
     computed per context, so the cache must be gated on non-alignable surfaces — which is why it
     isn't cached today. ODD here is non-aligned, so effectively all surfaces would benefit.

## Second pass

After the first pass (gen3 at ~4.10 s / 39.19e9 instructions vs gen1 ~4.02 s / 37.60e9), a fresh
`perf annotate` localized the remaining gen3-specific cost to four spots, each fixed and measured
independently:

1. **`AnyBase` zero-fill** — the default member initializer on the small buffer value-initialized
   128 bytes on every `pushState` (`rep stos`, ~29 % of `createState` self time). Removed; the
   buffer is only ever read through the handler. The default constructor became user-provided so
   `const` AnyBase objects can still be default-initialized.
2. **Per-step `isDefault()`** — paid a static-guard check plus a PLT call to
   `typeHash<EmptyState>()` every step (~7 % of `getNextTargetGen3` self time). The defaultness
   only changes at volume transitions, so the navigator now caches it as a bool in its state.
   Together with (1): 39.19 → 39.09e9 instructions, wall ~4.10 → ~4.06 s.
3. **Candidate copy-out** — `resolveCandidates` copied the path-length-filtered candidates into a
   separate `small_vector` per rebuild. The navigator now consumes the (sorted) stream directly
   and applies the near/far window lazily while advancing — identical acceptance in identical
   order, no copy, no spill, and never-consumed candidates skip the filter. Wall ~4.06 →
   ~3.94-4.00 s; instructions flat (the removed work was memory traffic).
4. **Dedup pre-pass + `checkPathLength`** — the O(n²) pre-intersection dedup only matters when
   external/free surfaces were appended, which the navigator knows; a defaulted
   `candidatesAreUnique` hint on `NavigationStream::initialize` skips it (post-sort unique remains
   as backstop, unit-tested). Independently, `detail::checkPathLength` was an out-of-line call
   doing two comparisons plus up to four verbose-filter checks, hot in **both** gens; the fast
   path is now inline with the verbose diagnostics out of line. Together: gen3 39.10 → 38.09e9
   (−2.6 %), **gen1 37.60 → 36.94e9 (−1.7 %)** — the first shared win in this series.
5. **Stateless-policy state skip** — the remaining `createState`/`popState` round-trip per volume
   entry (virtual child dispatch + `EmptyState` pushes/pops) is pure overhead when no policy
   carries state. Each `INavigationPolicy` now caches whether it pushes only default states,
   probed once at the end of `Blueprint::construct` under a documented contract (state
   *defaultness* must not depend on context/arguments — the type is part of the policy's
   identity). The property lives on the policy, not the volume, so `TrackingVolume` gains no
   member or API, and a policy attached later conservatively reports stateful until re-probed.
   The navigator skips create and pop under the same immutable condition; debug builds always
   exercise the subsystem and assert the contract. gen3 38.09 → 37.13e9 (−2.5 %), wall ~3.82 s,
   15.2 → 14.2 malloc/prop.

**Net after both passes** (100k ODD propagations, single thread):

| | malloc/prop | instructions | wall |
|------|------------|--------------|------|
| gen1 (unchanged except 4.) | 13.2 | 36.93e9 | ~4.0 s |
| gen3 before any work | 99.3 | — | ~4.55 s |
| gen3 after pass 1 | 15.2 | 39.19e9 | ~4.10 s |
| gen3 after pass 2 | **14.2** | **37.13e9** | **~3.82 s** |

Gen3 went from ~14 % slower than gen1 to **faster in wall clock** and within 0.5 % in
instructions. The remaining structural difference (gen3 re-intersects the volume candidate set on
volume entry, dominated by the virtual `Surface::intersect` dispatch inside
`NavigationStream::initialize`) is shared-shape work that gen1 pays in `compatibleSurfaces`
instead.

## Method (for reproduction)

A small driver runs the `PropagationAlgorithm` over a muon particle gun on the ODD built in
either mode (`getOpenDataDetector(gen3=…)`), with the same `EigenStepper`, constant 2 T field,
single thread, and IO disabled; the only difference between runs is the geometry + navigator.
Three independent signals were used:

- **Wall time** — mean of repeated runs at 100k propagations (100 events × 1000 tracks).
- **Instructions** — `perf stat -e instructions:u`, a load-independent, sub-0.1 %-stable metric
  used for the small (<3 %) call-overhead deltas that wall time can't resolve.
- **Allocations per propagation** — an `LD_PRELOAD` `malloc`/`free` counter, measured at two
  workloads; the difference isolates the per-propagation rate from fixed startup/geometry cost.
  A backtrace-attributing variant identified the specific allocating call sites.

Hot spots were localized with `perf record -F 4000` (flat self-time) and `--call-graph dwarf`
(caller trees). The measurement harness is not part of this change.

### Caveats
- Material is off; enabling it adds shared stepper-side cost to both modes and shifts the
  relative navigation fraction down, but does not change the gen3-vs-gen1 allocation delta.
- Single detector (ODD), single stepper (Eigen), constant field. Findings are about geometry
  navigation, which is the only thing that differs between the two modes here.
