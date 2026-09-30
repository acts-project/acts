@defgroup python_custom_algorithms Custom Algorithms
@ingroup python_bindings
@brief Writing your own `IAlgorithm` or `IReader` in Python, and moving data through the whiteboard.

Every step in a `Sequencer` — reader, algorithm, or writer — implements a small C++ interface.
Python subclasses of `acts.examples.IAlgorithm` and `acts.examples.IReader` are accepted directly
by `Sequencer.addAlgorithm()`/`addReader()`, alongside the built-in bindings, so you can write a
custom step without touching C++.

# A minimal algorithm

@snippet{trimleft} examples/test_custom_algorithm.py Custom algorithm filtering particles by pT

Points to note:

- The constructor must call the base class `__init__` (`acts.examples.IAlgorithm.__init__(self,
  name, level)`) before doing anything else.
- `execute(self, context)` is the only method you must override; it returns an
  `acts.examples.ProcessCode` (`SUCCESS`, `SKIP`, `ABORT`, or `END`).
- `self.logger` is available once the base constructor has run.

# Reading and writing the whiteboard

Data is passed between steps via an in-memory `WhiteBoard`, keyed by string. `ReadDataHandle`
and `WriteDataHandle` are the typed accessors:

```python
self.inputParticles = acts.examples.ReadDataHandle(
    self, acts.examples.SimParticleContainer, "InputParticles"
)
self.inputParticles.initialize("particles_generated")  # the whiteboard key to bind to
```

```python
self.outputParticles = acts.examples.WriteDataHandle(
    self, acts.examples.SimParticleContainer, "OutputParticles"
)
self.outputParticles.initialize("particles_high_pt")
```

Inside `execute`, read with `handle(context.eventStore)`, write with `handle(context, value)`.
`context.eventStore` also gives you `.exists(key)` and `.keys` directly, if you need to probe the
whiteboard rather than go through a typed handle.

> [!note]
> A `ReadDataHandle`/`WriteDataHandle` can only be constructed for a type that has been
> registered for whiteboard access on the C++ side (`SimParticleContainer`, `SimHitContainer`,
> `MeasurementContainer`, `ProtoTrackContainer`, `ClusterContainer`, `ConstTrackContainer`, space
> points, seeds, and a handful of others). Constructing one for an unregistered type raises
> `TypeError: ... is not registered for WhiteBoard access` immediately, not at first use.

# A custom reader

`IReader` follows the same pattern, but implements `availableEvents()` (returning
`(firstEvent, lastEventExclusive)`) instead of taking whiteboard input. `acts.examples.uproot`
ships two ready-to-use examples that read ROOT files without the ROOT plugin — read them for a
complete, working `IReader` implementation, including buffered reads across events.

# Building a `TrackContainer` by hand

Producing tracks in a custom algorithm (rather than reading them) means building a
`acts.examples.TrackContainer` directly:

```python
container = acts.examples.TrackContainer()
track = container.makeTrack()
track.parameters = acts.BoundVector(loc0, loc1, phi, theta, qOverP, time)
track.particleHypothesis = acts.ParticleHypothesis.muon  # defaults to pion, see caveats

for sourceLink in sourceLinksForThisTrack:
    trackState = track.appendTrackState()
    trackState.typeFlags.isMeasurement = True
    trackState.uncalibratedSourceLink = sourceLink
    trackState.referenceSurface = surfaceForSourceLink(sourceLink)

track.nMeasurements = len(sourceLinksForThisTrack)
self.outputTracks(context, container.makeConst())
```

`Examples/Scripts/Python/pypi_finding_fitting_demo.py` is the canonical worked example: a
complete, ROOT-free chain with a custom Python track finder (turning space points into
`ProtoTrack`s) and a custom Python track fitter (turning `ProtoTrack`s into a `TrackContainer`),
followed by truth matching and the ROOT-free performance writers from
@ref python_performance_plotting. It only uses what the PyPI wheel provides — copy it as a
starting point.

# Keeping the algorithm thin

For anything beyond a few lines, keep the `IAlgorithm`/`IReader` subclass itself minimal — base
class calls, data handles, and a call into your own plain Python class or function — rather than
putting the actual logic inline. This keeps the ACTS-specific plumbing (handles, `ProcessCode`,
context) separate from logic you likely want to unit test on its own.
