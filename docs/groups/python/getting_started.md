@defgroup python_getting_started Getting Started
@ingroup python_bindings
@brief The `Sequencer`, and the high-level `add*` helper API most users actually write against.

# The Sequencer

Every ACTS Python workflow is built around `acts.examples.Sequencer`: readers, algorithms, and
writers are added to it, and `s.run()` drives them event by event through an in-memory
`WhiteBoard` (see @ref python_custom_algorithms for how to read/write it directly).

Here is a minimal example that sets up particle propagation and runs a few events:

@snippet{trimleft} examples/test_generic.py Basic propagation example with GenericDetector

# The `add*` helper API

Wiring up algorithms one `Config` at a time is rarely how ACTS Python is actually written.
Instead, `acts.examples.simulation` and `acts.examples.reconstruction` provide chain-building
functions that each add a whole pipeline stage — generation, simulation, digitization, seeding,
track finding, fitting, vertexing — to a `Sequencer`, wiring the whiteboard keys between stages
for you. Nearly every script under `Examples/Scripts/Python/` is built by composing these.

```python
from pathlib import Path
import acts
from acts import UnitConstants as u
import acts.examples
from acts.examples.simulation import addParticleGun, addFatras, addDigitization, MomentumConfig, EtaConfig, ParticleConfig
from acts.examples.reconstruction import addSeeding, addKalmanTracks

s = acts.examples.Sequencer(events=100, numThreads=-1)
rnd = acts.examples.RandomNumbers(seed=42)

addParticleGun(
    s,
    MomentumConfig(1 * u.GeV, 10 * u.GeV, transverse=True),
    EtaConfig(-2.0, 2.0),
    ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
    rnd=rnd,
)
addFatras(s, trackingGeometry, field, rnd=rnd)
addDigitization(s, trackingGeometry, field, digiConfigFile=digiConfigFile, rnd=rnd)
addSeeding(s, trackingGeometry, field, geoSelectionConfigFile=geoSelectionConfigFile)
addKalmanTracks(s, trackingGeometry, field)

s.run()
```

Each function follows the same shape: it takes the `Sequencer` plus whatever inputs it needs
(geometry, field, config files, a shared `RandomNumbers`), reads its inputs from the whiteboard
keys the previous stage wrote (with sensible defaults you can override), and returns after
appending its algorithms — nothing runs until `s.run()`.

> [!tip]
> These functions are the best reference for the underlying algorithm `Config` objects: reading
> `Python/Examples/python/simulation.py` and `reconstruction.py` shows exactly which config
> fields matter and which whiteboard keys are conventionally used, even if you end up building
> a chain by hand instead of calling them.

Full working chains built this way: `Examples/Scripts/Python/full_chain_odd.py` (the reference
OpenDataDetector chain), `truth_tracking_kalman.py` (truth-seeded, the simplest end-to-end
example), and `ckf_tracks.py`.

# Config objects and keyword arguments

Every bound `Config` class accepts its fields as keyword arguments directly (`Algorithm(level=...,
inputFoo="x", outputBar="y")` instead of building a `Config` object by hand first). This comes
from a Python-side patch applied to every binding module (`Python/Core/python/_adapter.py`), not
from pybind11 itself. See @ref python_caveats for how this interacts with nested config fields.
