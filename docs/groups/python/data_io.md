@defgroup python_data_io Data Reading and Writing
@ingroup python_bindings
@brief The reader/writer inventory, including the ROOT-free path used by the PyPI wheel.

# Reader/writer inventory

Every reader and writer is a `SequenceElement` added via `s.addReader(...)`/`s.addWriter(...)`.
Formats available depend on which plugins were built (see @ref python_installation for what the
PyPI wheel includes):

- **CSV** (`acts.examples`) — always available: `Csv{Particle,Measurement,SimHit,SpacePoint,
  Track,TrackParameter,ProtoTrack}{Reader,Writer}`, plus a few detector-specific ones.
- **ROOT** (`acts.examples.root`, source builds only) — `Root{Particle,Vertex,SimHit,
  TrackSummary,MaterialTrack,...}{Reader,Writer}`.
- **EDM4hep/podio** (`acts.examples.edm4hep`, source builds only) — `PodioReader`/`PodioWriter`.
- **HepMC3** (`acts.examples.hepmc3`) — `HepMC3Reader`/`HepMC3Writer`.
- **Arrow/Parquet** (`acts.examples.arrow`) — **built into the PyPI wheel**, see below.
- **JSON** (`acts.json`, `acts.examples.json`) — geometry, material maps, and lookup tables, not
  event data.

# The ROOT-free path (PyPI wheel)

## Arrow/Parquet

`acts.examples.arrow.ParquetReader`/`ParquetWriter` read and write whole event collections as
sharded Parquet files, with an explicit expected schema:

```python
import acts.arrow  # module-level table schemas, e.g. acts.arrow.particleSchema()
import acts.examples.arrow

reader = acts.examples.arrow.ParquetReader(
    level=acts.logging.INFO,
    inputDir=str(inputDir),
    collections={"particles_arrow": "particles_generated"},
    expectedSchemas={"particles_arrow": acts.arrow.particleSchema()},
)
s.addReader(reader)
```

`ColliderMLRelease1InputConverter` reads the [ColliderML](https://huggingface.co/CERN) Release 1
Parquet schema directly into `particles`, `simhits`, `measurements`, and the associated index
maps, including a geometry-ID remapping CSV (`geoIdMapPath`) between the dataset's geometry and
your `TrackingGeometry`.

## Uproot readers

`acts.examples.uproot` provides `UprootParticleReader` and `UprootSimHitReader`: pure-Python
`IReader`s that read the exact file format written by `RootParticleWriter`/`RootSimHitWriter`,
using [uproot](https://uproot.readthedocs.io/) instead of the ROOT plugin — useful when you have
ROOT files from elsewhere but don't have (or want) a ROOT-enabled ACTS build. They double as a
complete, real-world example of a custom `IReader` (see @ref python_custom_algorithms), including
buffered multi-event reads.

## ROOT-free fitting and evaluation

`acts.examples.scipy.makeScipyHistogramFitFunction()` provides a Gaussian-fit backend for the
performance writers below, replacing the ROOT-based fit used by `ActsPlugins::RootHistogramFit`.
See @ref python_performance_plotting.

# Geometry without DD4hep

Load a `TrackingGeometry` from a JSON dump instead of building it from DD4hep or Geant4 — this is
how you get a real detector geometry into a PyPI-only script:

```python
import acts.json

gctx = acts.GeometryContext()
converter = acts.json.TrackingGeometryJsonConverter()
trackingGeometry = converter.fromFile(gctx, "geometry.json")
```

Material maps can be attached the same way with `acts.json.JsonMaterialDecorator`.
`Examples/Scripts/Python/geometry.py` shows how to produce such a JSON dump from an existing
`TrackingGeometry` (`toFile`/`toJson`), for example from a source build with DD4hep enabled.
