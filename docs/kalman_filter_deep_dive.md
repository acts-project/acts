# The Kalman Filter in ACTS — A Deep Dive

> A walk-through of how the Kalman track fitter is implemented in this
> repository and how the pieces fit together. File/line references point at
> the sources as of the `claude/kalman-filter-deep-dive` branch.

## 1. What ACTS calls "the Kalman filter"

In ACTS the Kalman filter shows up in three related-but-distinct places. It
helps to keep them separate:

| Component | Purpose | Location |
|---|---|---|
| **`KalmanFitter`** | Track *fitting*: hits already assigned to a track, estimate the parameters | `Core/include/Acts/TrackFitting/KalmanFitter.hpp` |
| **`CombinatorialKalmanFilter` (CKF)** | Track *finding*: simultaneously associates hits and fits, branching on ambiguities | `Core/include/Acts/TrackFinding/CombinatorialKalmanFilter.hpp` |
| **`GaussianSumFitter` (GSF)** | A multi-component extension of the KF for non-Gaussian energy loss (electrons) | `Core/include/Acts/TrackFitting/GaussianSumFitter.hpp` |

This document focuses on the **`KalmanFitter`**, the canonical implementation,
and the shared math kernels (`GainMatrixUpdater`, `GainMatrixSmoother`,
`MbfSmoother`) that the CKF and GSF reuse.

The single most important architectural idea: **the Kalman filter is not a
standalone loop. It is a *plugin* (an "Actor") that rides along inside the
track propagator.** The propagator transports the track through the detector
geometry surface-by-surface; every time it lands on a surface, it calls back
into the Kalman actor, which does the measurement update. This is the key to
understanding everything else.

## 2. The big picture: the KF as a propagator Actor

ACTS's `Propagator` walks a track through the geometry using a *stepper*
(numerical integration of the equations of motion in the B-field, e.g.
`SympyStepper`) and a *navigator* (decides which surface comes next). The
propagator supports a list of **Actors** — callbacks invoked at each step. The
`KalmanFitter` defines a private `Actor` class (`KalmanFitter.hpp:286`) that is
injected into the propagator.

```text
Propagator.propagate()
   └─ for each step:
        stepper integrates ODE to next surface
        navigator identifies currentSurface
        Actor::act(state, stepper, navigator, result)   ← the Kalman logic lives here
        Actor::checkAbort(...)                           ← stop when track complete / target reached
```

So a "fit" is really a *propagation with a measurement-updating actor
attached*. The fitter object itself (`KalmanFitter`, `KalmanFitter.hpp:252`) is
thin — it owns a `propagator_t m_propagator` and two loggers
(`KalmanFitter.hpp:273-278`), and its real job is to set up propagator options,
run the propagation, then orchestrate smoothing.

The class is templated on `<propagator_t, traj_t>`:

- `propagator_t` — the propagator type (which fixes the stepper and navigator).
- `traj_t` — the trajectory storage backend, in practice
  `Acts::VectorMultiTrajectory`.

A `static constexpr bool isDirectNavigator` (`KalmanFitter.hpp:258`) selects
between the standard `Navigator` (figures out surfaces from geometry on the
fly) and the `DirectNavigator` (you hand it the exact surface sequence — used
for refitting).

## 3. The data structures that flow through a fit

### `KalmanFitterOptions<traj_t>` (`KalmanFitter.hpp:108-186`)

Everything that steers one fit:

- The three **contexts**: `geoContext`, `magFieldContext`,
  `calibrationContext` — ACTS's mechanism for thread-safe, per-event
  conditions data (alignment, field, calibration).
- `extensions` — the customization delegates (see §9).
- `propagatorPlainOptions` — step sizes, tolerances, path limits.
- `referenceSurface` + `referenceSurfaceStrategy` — where you want the final
  parameters reported (e.g. the perigee/beamline).
- Physics toggles: `multipleScattering`, `energyLoss`, `freeToBoundCorrection`.
- Smoothing controls: `reverseFiltering`, `reverseFilteringCovarianceScaling`
  (default `100.0`).

### `KalmanFitterExtensions<traj_t>` (`KalmanFitter.hpp:43-103`)

A bundle of `Acts::Delegate`s (type-erased function pointers, no heap
allocation) that let you swap behavior without recompiling the fitter:

- `calibrator` — turns a raw `SourceLink` into a calibrated measurement on the
  track state.
- `updater` — *the Kalman update itself* (default: `GainMatrixUpdater`).
- `smoother` — backward smoothing (default in examples: `MbfSmoother`).
- `outlierFinder` — decides whether a measurement is an outlier (skip the
  update, mark the state).
- `reverseFilteringLogic` — decides per-track whether to smooth via reverse
  filtering vs. the smoother.
- `surfaceAccessor` — maps a `SourceLink` to its `Surface`.

The default constructor wires every delegate to a `void` stub
(`KalmanFitter.hpp:94-102`) — these throw or no-op, and exist so the type is
constructible for unit tests. **You must connect real implementations before
fitting** (the examples layer does this, §10).

### `KalmanFitterResult<traj_t>` (`KalmanFitter.hpp:189-229`)

The mutable scratchpad the actor writes into during propagation:

- `fittedStates` — pointer to the `MultiTrajectory` where track states are
  appended.
- `lastMeasurementIndex` / `lastTrackIndex` — tips of the trajectory linked
  list.
- Counters: `measurementStates`, `measurementHoles`, `processedStates`.
- `missedActiveSurfaces` — sensitive surfaces crossed with no hit (holes).
- `fittedParameters` — parameters bound to the target surface once reached.
- `finished` — flips true to terminate propagation.

### The `MultiTrajectory` and track states

Track states are stored in a `MultiTrajectory` (a columnar, index-linked
container). Each state can hold several parameter vectors, selected by a
**`TrackStatePropMask`** bitmask — `Predicted`, `Filtered`, `Smoothed`,
`Jacobian`, `Calibrated`. These three parameter sets are the heart of the
Kalman formalism:

- **predicted** — extrapolation onto this surface using only *past*
  measurements.
- **filtered** — predicted, updated with *this* surface's measurement (a
  weighted mean).
- **smoothed** — the best estimate using *all* measurements, computed in a
  backward pass.

States also carry **type flags** (`setHasParameters`, `setHasMaterial`,
`setIsMeasurement`, `setIsOutlier`, `setIsHole`) that downstream code (and
`calculateTrackQuantities`) reads.

## 4. The `fit()` entry point and the two-pass structure

Public `fit(...)` overloads (`KalmanFitter.hpp:667` and `:695`) just forward to
`fit_impl` (`KalmanFitter.hpp:835`). The flow:

```text
fit_impl
 ├─ make_propagator_options(... reverseDirection=false)   build forward options
 ├─ filter_impl(...)                                       FORWARD FILTER PASS
 │     └─ propagator.propagate()  →  Actor runs at every surface
 ├─ seed last measurement: smoothed ← filtered
 ├─ decide: reverse-filter or smooth?
 │     ├─ reverse filtering:  filter_impl(... reverseDirection=true)   SECOND PROPAGATION
 │     └─ direct smoothing:   extensions.smoother(...)                 backward sweep, no propagation
 ├─ extrapolate to referenceSurface (if not already there)
 └─ calculateTrackQuantities(track)
```

Two key setup details inside `make_propagator_options`
(`KalmanFitter.hpp:710`):

1. **Measurements become a map.** Input source links are copied into
   `std::unordered_map<const Surface*, SourceLink> inputMeasurements` keyed by
   surface (`KalmanFitter.hpp:726-731`). The fit therefore **does not require
   measurements to be sorted** — ordering emerges from navigation. This is
   stated explicitly in the class comment (`KalmanFitter.hpp:247-249`).
2. **Surfaces are registered with the navigator.** For the standard navigator,
   each measurement surface is appended as an *external surface*
   (`appendExternalSurface`, `KalmanFitter.hpp:744-749`) so the navigator
   targets them even if a boundary check would otherwise miss. For the
   `DirectNavigator`, the full surface sequence is set directly
   (`KalmanFitter.hpp:750-755`).

The actor is then populated with the measurement map, target surface, physics
toggles, extensions, and logger (`KalmanFitter.hpp:758-766`).

`filter_impl` (`KalmanFitter.hpp:771`) makes a propagator state, points
`kalmanResult.fittedStates` at the track-state container, runs `propagate`,
and — crucially — fails the fit if **zero** measurement states were produced
(`KalmanFitter.hpp:801-804`, `KalmanFitterError::NoMeasurementFound`). On
success it materializes a `TrackProxy`, sets its tip index, copies fitted
parameters/covariance, and calls `calculateTrackQuantities`.

## 5. The forward filter pass — `Actor::act` → `Actor::filter`

### `act` (`KalmanFitter.hpp:339`)

Called at every propagation step. It:

1. Returns early if `result.finished`.
2. On the first call, arms the loop-protection path limit
   (`setupLoopProtection`, `KalmanFitter.hpp:355-359`) — prevents a
   low-momentum track from looping forever in the field.
3. If the navigator reports a `currentSurface`, calls `filter(...)` on it
   (`KalmanFitter.hpp:363-379`).
4. Checks the battery of termination conditions
   (`KalmanFitter.hpp:383-394`): track complete
   (`measurementStates == inputMeasurements.size()`), end of world, volume
   constraint, path limit, target reached. If the target surface is reached, it
   binds the final parameters there (`stepper.boundState`,
   `KalmanFitter.hpp:407-415`) and sets `finished = true`.

`checkAbort` (`KalmanFitter.hpp:426`) just returns `result.finished` — that's
how the propagator knows to stop.

### `filter` (`KalmanFitter.hpp:445`) — the core per-surface logic

This is where the actual Kalman bookkeeping happens. It branches into three
cases based on what the surface *is*:

**Case A — the surface carries a measurement**
(`inputMeasurements.find(&surface)` hits, `KalmanFitter.hpp:453`):

1. **Transport covariance to the surface**:
   `stepper.transportCovarianceToBound` (`KalmanFitter.hpp:459`).
2. **Pre-update material effects**:
   `performMaterialInteraction(... PreUpdate ...)`
   (`KalmanFitter.hpp:463-468`). This inflates the covariance with
   multiple-scattering noise and shifts q/p for energy loss, *before* the
   measurement update. (`determineMaterialUpdateMode` avoids double-counting
   material at start/target surfaces.)
3. **Create a track state** with mask
   `Predicted | Filtered | Jacobian | Calibrated` (`KalmanFitter.hpp:471-475`).
4. **Bind the propagated state**: `stepper.boundState` gives
   `(boundParams, jacobian, pathLength)`. These fill `predicted()`,
   `predictedCovariance()` (taken from `state.stepping.cov`), `jacobian()`,
   `pathLength()` (`KalmanFitter.hpp:494-498`). **The Jacobian stored on a state
   is the transport Jacobian from the *previous* surface to *this* one** — the
   smoother relies on this convention.
5. **Calibrate**: `extensions.calibrator(...)` turns the source link into a
   calibrated measurement now that predicted parameters exist
   (`KalmanFitter.hpp:502`).
6. **Outlier test**: `extensions.outlierFinder(...)`
   (`KalmanFitter.hpp:518`).
   - **Not an outlier** → run `extensions.updater(...)` (the Kalman update,
     §6), tag `setIsMeasurement` (`KalmanFitter.hpp:520-527`).
   - **Outlier** → don't update; tag `setIsOutlier`; alias `Filtered` to
     `Predicted` via `shareFrom` so the filtered slot is still readable
     (`KalmanFitter.hpp:528-537`).
7. **Update the stepper** with the filtered state (`stepper.update`,
   `KalmanFitter.hpp:547-551`) so subsequent propagation continues from the
   *improved* estimate — this is what makes it a recursive filter. Increment
   `measurementStates`.
8. **Post-update material effects**: a second
   `performMaterialInteraction(... PostUpdate ...)`
   (`KalmanFitter.hpp:557-562`). Splitting material pre/post the measurement is
   how ACTS models material on the same surface as the hit.

**Case B — no measurement, but a hole or material surface**
(`(precedingMeasurementExists && surfaceIsSensitive) || surfaceHasMaterial`,
`KalmanFitter.hpp:572`):

Creates a `Predicted | Jacobian` state, binds it, aliases
`Filtered ← Predicted` (`shareFrom`, `KalmanFitter.hpp:604-606`), and tags it:

- **Hole** (`setIsHole`) if it's a sensitive surface after the first
  measurement but with no hit → recorded in `missedActiveSurfaces`
  (`KalmanFitter.hpp:616-633`).
- Sensitive-but-before-first-measurement → skipped silently (no hole counted
  before the track starts).
- Pure material surface → `setHasMaterial`, used to apply scattering/energy-loss
  (`FullUpdate` mode, `KalmanFitter.hpp:638-643`).

**Case C** — neither (e.g. a passive surface before the first hit) → nothing
recorded.

This three-way split is why the resulting trajectory faithfully records
measurements, outliers, holes, and material — the information later quality
cuts and the smoother need.

## 6. The Kalman update math — `GainMatrixUpdater`

The default `updater` delegate. Header `GainMatrixUpdater.hpp`, dispatch in
`GainMatrixUpdater.cpp`, math in `detail/GainMatrixUpdaterImpl.hpp`.

`operator()` (`GainMatrixUpdater.hpp:38`) asserts the state has calibrated,
predicted, and (allocated) filtered slots, then calls `visitMeasurement`.
Because a measurement can be 1–6 dimensional (pixel = 2, strip = 1, etc.), the
dimension `N` is only known at runtime. `visitMeasurement`
(`GainMatrixUpdater.cpp:19`) uses `visit_measurement(calibratedSize(), ...)` to
dispatch to a compile-time-sized `visitMeasurementImpl<N>` — the templates are
explicitly instantiated for N = 1..6 (`GainMatrixUpdaterImpl.hpp:99-104`) to
keep compile times sane.

The math (`GainMatrixUpdaterImpl.hpp:21-91`), with $x^-, P^-$ = predicted,
$m, R$ = calibrated measurement + covariance, $H$ = projector from bound
(6-dim) to measurement (N-dim) space:

**Gain matrix** (`:52`):

$$ K = P^- H^\top \left( H P^- H^\top + R \right)^{-1} $$

**Filtered parameters** (`:65`, then phi/theta normalized):

$$ x^+ = x^- + K\left( m - H x^- \right) $$

**Filtered covariance** (`:69-75`) — standard form by default, optional
**Joseph form** (`m_useJosephFormulation`) for numerical stability:

$$ P^+ = (I - KH)\,P^- \qquad\text{or}\qquad P^+ = (I-KH)P^-(I-KH)^\top + K R K^\top $$

**Chi-square contribution** of this state (`:80-88`), using the *filtered*
residual:

$$ r = m - H x^+, \qquad \chi^2 = r^\top \big[(I - HK)R\big]^{-1} r $$

If $K$ contains NaNs the update fails with `KalmanFitterError::UpdateFailed`.
The projector $H$ comes from a `FixedBoundSubspaceHelper` built from the
state's stored subspace indices (`:36-43`) — i.e. which of the 6 bound
coordinates the detector actually measures.

## 7. Smoothing — three ways to get the "best" estimate

The forward filter gives `filtered` = best estimate at each surface using only
past+present hits. To incorporate *future* hits you smooth. The KF seeds the
last measurement's `smoothed ← filtered` (`KalmanFitter.hpp:865-867`) — at the
track end there is no future, so smoothed = filtered — and then propagates that
information backward. ACTS offers two backward-smoother kernels and one
alternative full-reverse-fit.

### 7a. `GainMatrixSmoother` — Rauch–Tung–Striebel (RTS)

File: `GainMatrixSmoother.hpp` + `GainMatrixSmoother.cpp`.

`operator()` (`GainMatrixSmoother.hpp:43`) seeds the entry state
(`Smoothed ← Filtered`) and walks backward via `trajectory.applyBackwards`
(`:69`), calling `calculate` for each state. In `calculate`
(`GainMatrixSmoother.cpp:18`), `ts` is the state being smoothed and `prev_ts`
is the *already-smoothed* neighbor toward the track end. Since
`prev_ts.jacobian()` is the forward transport from `ts` to `prev_ts`:

**Smoother gain** (`GainMatrixSmoother.cpp:32`):

$$ G = C^f_k\, F_k^\top \left(P^-_{k+1}\right)^{-1} $$

**Smoothed parameters** (`:65`) and **covariance** (`:75`):

$$ x^s_k = x^f_k + G\left(x^s_{k+1} - x^p_{k+1}\right) $$

$$ C^s_k = C^f_k + G\left(C^s_{k+1} - P^-_{k+1}\right) G^\top $$

It requires inverting the full $6\times6$ predicted covariance at each step. An
optional `doCovCheckAndAttemptFix` validates semi-positive-definiteness
(`:79-92`).

### 7b. `MbfSmoother` — Modified Bryson–Frazier

File: `MbfSmoother.hpp` + `.cpp` + `detail/MbfSmootherImpl.hpp`. **This is the
default in the examples** (`KalmanFitterFunction.cpp:77`).

The MBF smoother propagates two adjoint accumulators backward instead of
inverting the full covariance — it only needs the inverse of the (small, N×N)
residual covariance, which the filter could even cache. It tracks
$\hat\Lambda$ (`bigLambdaHat`, a $6\times6$ matrix) and $\hat\lambda$
(`smallLambdaHat`, a 6-vector), both starting at zero
(`MbfSmoother.hpp:61-62`).

For each state, backward (`MbfSmoother.hpp:64-86`):

**Smoothed estimate** (`MbfSmoother.cpp:16-27`):

$$ x^s = x^f - C^f\,\hat\lambda, \qquad C^s = C^f - C^f\,\hat\Lambda\,C^f $$

Then update the accumulators depending on state type:

- **Measurement state** (`MbfSmootherImpl.hpp:17-58`): with residual covariance
  $S = HP^-H^\top + R$, gain $K = P^-H^\top S^{-1}$, and $\hat C = I - KH$:

  $$ \tilde\Lambda = H^\top S^{-1} H + \hat C^\top \hat\Lambda \hat C, \qquad \tilde\lambda = -H^\top S^{-1} y + \hat C^\top \hat\lambda $$

  then transport across the Jacobian $F$: $\hat\Lambda = F^\top \tilde\Lambda F$,
  $\hat\lambda = F^\top\tilde\lambda$.
- **Non-measurement state** (hole/material, `MbfSmoother.cpp:29-36`): just
  transport — $\hat\Lambda = F^\top \hat\Lambda F$,
  $\hat\lambda = F^\top \hat\lambda$.

MBF is generally preferred: cheaper (no full-matrix inverse) and numerically
nicer. The class comment links the Wikipedia derivation
(`MbfSmoother.hpp:36-38`).

### 7c. Reverse filtering — a full second propagation

Chosen when `kfOptions.reverseFiltering` is set or `reverseFilteringLogic(...)`
returns true for this track (`KalmanFitter.hpp:871-874`). Used mostly for
low-momentum tracks where linearized smoothing is inaccurate.

Instead of an analytic backward sweep, it **runs the whole filter again in the
reverse direction** (`make_propagator_options(... reverseDirection=true)`,
`filter_impl`, `KalmanFitter.hpp:882-885`). The reverse pass starts from the
forward fit's last measurement, with its covariance **inflated by
`reverseFilteringCovarianceScaling`** (default 100×,
`KalmanFitter.hpp:880-881`) so the forward information doesn't dominate. After
the reverse fit, the first measurement's `smoothed` is taken from the reverse
pass's last-measurement `filtered` (`KalmanFitter.hpp:908-910`), a consistency
check guards against mismatched reference surfaces
(`KalmanFitter.hpp:900-907`), the final parameters are taken from the reverse
track, and the temporary reverse track is removed
(`KalmanFitter.hpp:912-919`).

The fit records which path it took in the `"smoothed"` / `"reversed"` boolean
columns if present (`KalmanFitter.hpp:948-955`).

## 8. Finishing up — reference surface and track quantities

If the chosen path didn't already leave parameters on the requested
`referenceSurface`, the fit does a final extrapolation there via
`extrapolateTrackToReferenceSurface` (`KalmanFitter.hpp:933-946`), using
`referenceSurfaceStrategy` (e.g. `first` for the perigee).

`calculateTrackQuantities` (`TrackHelpers.hpp:388`) then aggregates the
per-state results into track-level summaries by sweeping the states in
reverse: counts holes/outliers/measurements/shared hits, and sums **only
measurement states'** chi² into `track.chi2()` and their calibrated sizes into
`track.nDoF()` (`TrackHelpers.hpp:399-414`).

## 9. The extension points (and why they're delegates)

The fitter's flexibility comes from the six delegates in
`KalmanFitterExtensions` (§3). They're `Acts::Delegate`s — type-erased,
allocation-free, connected at runtime. Default "void" stubs live in
`detail/VoidFitterComponents.hpp`:

- `voidFitterUpdater` (`:29-36`) — copies predicted→filtered (no measurement
  update); used to test pure propagation.
- `voidFitterSmoother` (`:38-48`) — copies filtered→smoothed backward.
- `voidOutlierFinder` (`:50-53`) — never an outlier.
- `voidReverseFilteringLogic` (`:55-59`) — never reverse-filter.
- `voidFitterCalibrator` / `voidSurfaceAccessor` (`:21-27`, `:61-63`) —
  **throw** if called (you *must* supply real ones).

This delegate design is what lets the *same* `KalmanFitter` code drive a pixel
detector, a straw-tube tracker, an alignment refit, or a unit test — only the
connected functions change.

## 10. How it's actually wired up (the Examples layer)

`Examples/Algorithms/TrackFitting/src/KalmanFitterFunction.cpp` shows
production wiring:

- Stepper = `SympyStepper`; two fitters built — one with the standard
  `Navigator`, one with `DirectNavigator` for refits (`:39-44`, `:170-190`).
- `makeKfOptions` (`:94`) connects the real delegates: `GainMatrixUpdater` as
  `updater`, **`MbfSmoother` as `smoother`** (`:98-103`), a
  `SimpleReverseFilteringLogic` that triggers reverse filtering below a momentum
  threshold (`:50-58`, `:104-106`), and a `SimpleOutlierFinder` that cuts on
  predicted χ² (`:60-68`, `:107-108`).
- The calibrator and surface accessor are connected per-call, switching to a
  `RefittingCalibrator` when `doRefit` is set (`:119-131`).
- `useJosephFormulation` flows into the updater constructor (`:203`).

So the public-facing "Kalman fitter" a user runs is: `SympyStepper` +
`Navigator` + `KalmanFitter` actor + `GainMatrixUpdater` + `MbfSmoother`, with
chi²-based outlier rejection and momentum-gated reverse filtering.

## 11. How the pieces fit together (mental model)

```text
                    KalmanFitter.fit()
                          │
        ┌─────────────────┴───────────────────┐
        │  FORWARD PASS (propagation + actor)  │
        │                                      │
   Propagator ── Stepper (B-field ODE)         │
        │   └── Navigator (next surface)       │
        │                                      │
   Actor::act ──► Actor::filter(surface):      │
        │   1. transport cov to surface        │
        │   2. material PRE-update  ───────────┼── PointwiseMaterialInteraction
        │   3. predicted = boundState          │     (scattering noise, energy loss)
        │   4. calibrate  ─────────────────────┼── extensions.calibrator
        │   5. outlier?  ──────────────────────┼── extensions.outlierFinder
        │   6. UPDATE: filtered ───────────────┼── extensions.updater = GainMatrixUpdater
        │   7. stepper.update(filtered)        │       K, x⁺, P⁺, χ²
        │   8. material POST-update            │
        │   (holes / material handled too)     │
        └──────────────────┬───────────────────┘
                           │ seed: last meas. smoothed ← filtered
        ┌──────────────────┴───────────────────┐
        │  SMOOTHING (pick one)                 │
        │   • MbfSmoother (default)  ───────────┼── adjoint Λ̂, λ̂ backward sweep
        │   • GainMatrixSmoother (RTS)          │
        │   • Reverse filtering  ───────────────┼── second full propagation, inverted dir
        └──────────────────┬───────────────────┘
                           │
              extrapolate to referenceSurface
                           │
              calculateTrackQuantities (χ², nDoF, holes…)
                           │
                       TrackProxy
```

**Key takeaways:**

1. The KF is a **propagator actor**, not a standalone loop — geometry
   navigation drives surface ordering, so measurements need not be pre-sorted.
2. Every surface produces a **track state** carrying `predicted` / `filtered` /
   `smoothed` parameters; type flags distinguish measurement / outlier / hole /
   material.
3. The **gain-matrix update** (`GainMatrixUpdater`) is the recursive forward
   step; the stepper is re-seeded with the filtered state so the filter is
   genuinely recursive.
4. **Material** is modeled inline via pre/post-update interactions on each
   surface.
5. **Smoothing** runs the future information backward — cheaply via MBF/RTS, or
   robustly via a full reverse filter for low-momentum tracks.
6. Everything pluggable is a **`Delegate`** — the same code serves fitting,
   refitting, and (with branching) the CKF and (with mixtures) the GSF.
