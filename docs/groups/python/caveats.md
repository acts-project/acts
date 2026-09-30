@defgroup python_caveats Best Practices and Caveats
@ingroup python_bindings
@brief By-design behavior that is easy to get wrong. Binding gaps and bugs are tracked as
issues instead — see the note at the end of this page.

# Threading and the GIL

`Sequencer.run()` releases the Python GIL for the duration of the run
(`Python/Examples/src/Framework.cpp`, `run()`); every override point on a Python
`IAlgorithm`/`IReader`/`SequenceElement` (`execute`, `read`, `initialize`, `finalize`, `name`)
re-acquires it individually before calling back into Python. This lets C++-only algorithms run
fully in parallel across `numThreads`, but any *Python* algorithm effectively serializes on the
GIL, and if it holds non-thread-safe state (e.g. a CUDA model, a non-reentrant library handle),
run the `Sequencer` with `numThreads=1`.

# Config objects: nested fields don't mutate in place

Every bound `Config` field is exposed via `def_readwrite`
(`Python/Utilities/include/ActsPython/Utilities/Macros.hpp:45-46`). For a nested config struct,
reading `writer.config.resPlotToolConfig` returns a **copy**; mutating that copy has no effect
on `writer`. Read, modify, and reassign:

```python
res_plot_cfg = cfg.resPlotToolConfig
res_plot_cfg.varBinning["loc0"] = acts.Axis.regular(50, -1 * u.mm, 1 * u.mm, "d0 [mm]")
cfg.resPlotToolConfig = res_plot_cfg  # must reassign — in-place edits are silently discarded
```

# `particleHypothesis` defaults to pion

A newly created track (`TrackContainer.makeTrack()`) starts with `particleHypothesis` set to
pion (`Core/src/EventData/VectorTrackContainer.cpp:48`), regardless of what you're actually
fitting. If you build a `TrackContainer` by hand (@ref python_custom_algorithms), set
`track.particleHypothesis` explicitly.

# `TrackTruthMatcher` silently skips states without a real source link

Truth matching walks a track's states and skips any state that either isn't flagged
`isMeasurement`, is flagged `isOutlier`, or whose `uncalibratedSourceLink` doesn't hold an
`IndexSourceLink` (`Examples/Framework/src/Validation/TrackClassification.cpp:102-114`) —
silently, with no diagnostic. A custom track finder/fitter (@ref python_custom_algorithms) must
set a real `uncalibratedSourceLink` on every kept track state, or that state is invisible to
truth matching and performance evaluation.

# Performance writers scale with the *input* particle collection, not the matched tracks

`PatternRecognitionPerformanceCollector::fill()` loops over the whole `inputParticles` collection
you pass it once per track for a nested search
(`Examples/Framework/src/Validation/PatternRecognitionPerformanceCollector.cpp:152,177`). Pass an
already-selected/filtered particle collection (e.g. the output of `addDigiParticleSelection`),
not the raw generated one — passing the unfiltered collection can turn a few-second writer step
into minutes, with no error, just a slow run.

# Digitization silently drops hits that don't land on their surface

Smearing digitization looks up local coordinates on the hit's surface and, if that fails,
logs at `DEBUG` (invisible at the default `INFO` level) and skips the hit entirely
(`Examples/Algorithms/Digitization/src/DigitizationAlgorithm.cpp:293-299`), counting it in an
internal `skippedHits` tally that isn't surfaced anywhere by default. This can happen for
surfaces whose placement doesn't exactly match the geometry the hits were simulated against
(e.g. after a lossy geometry round-trip), and the position check uses a fixed absolute tolerance
(`Acts::s_onSurfaceTolerance = 1e-4`, `Core/include/Acts/Definitions/Tolerance.hpp:23`), not a
relative one. If your reconstructed efficiency is unexpectedly low, raise the log level to
`DEBUG` and check for this message before looking anywhere else.

# Splicing whiteboard keys with `addWhiteboardAlias`

The `add*` helpers read and write whiteboard keys with fixed conventions
(e.g. `addKalmanTracks` defaults to reading `truth_particle_tracks`). To feed a helper from a
step that used a different name, use `s.addWhiteboardAlias(newKey, existingKey)`
(`Python/Examples/src/Framework.cpp:383`) rather than renaming things upstream — see
`Examples/Scripts/Python/truth_tracking_kalman.py` for real usage.

# Renamed classes

`PythonTrackFinderPerformanceWriter` and `PythonTrackFitterPerformanceWriter` are deprecated
aliases (`Python/Examples/python/__init__.py:18-44`) for
`PythonPatternRecognitionPerformanceWriter` and `PythonTrackParameterPerformanceWriter`. The same
renaming applies to the ROOT-backed writers. Prefer the new names; the old ones still work but
raise a `DeprecationWarning`.

# Optional plugins fail loudly, on purpose

Importing a plugin submodule that wasn't built (e.g. `acts.dd4hep` without
`ACTS_BUILD_PLUGIN_DD4HEP=ON`) raises `ModuleNotFoundError`, with a warning naming the CMake flag
to build it (`Python/Plugins/python/plugin.py.in`). This is deliberate — treat it as "not built",
not as a bug to work around.

---

> [!note]
> A few behaviors reported against the Python bindings are binding **gaps**, not intentional
> design, and are tracked as issues rather than documented as permanent: `Surface::globalToLocal`
> is not exposed to Python (only `localToGlobal`); `BoundMatrix` cannot be filled from Python
> (only `.Identity()`/`.Zero()`; `np.asarray()` on it silently produces a useless 0-d object
> array); `ConstTrackStateProxy` (what you get reading an already-const track container) does
> not expose `uncalibratedSourceLink`/`referenceSurface`; and a container written through a
> `WriteDataHandle` is disowned on the Python side, so it must be read back through a
> `ReadDataHandle` in a later step rather than kept as a live Python reference.
