@defgroup python_data_io Data Reading and Writing
@ingroup python_bindings
@brief Readers and writers available in the PyPI wheel and optional source builds.

# Reader/writer inventory

Add readers with `s.addReader(...)` and writers with `s.addWriter(...)`. Availability depends on
the build configuration (see @ref python_bindings):

- **Only in a full installation, when enabled:**
  - ROOT (`acts.examples.root`): readers and writers for particles, sim hits, tracks, vertices,
    material, and performance output.
  - EDM4hep/podio (`acts.examples.edm4hep`): `PodioReader` and `PodioWriter`.
- **Also in the PyPI wheel:**
  - CSV (`acts.examples`): readers and writers for common event collections.
  - Arrow/Parquet (`acts.examples.arrow`): `ParquetReader` and `ParquetWriter`.
  - HepMC3 (`acts.examples.hepmc3`): `HepMC3Reader` and `HepMC3Writer`.
  - Uproot (`acts.examples.uproot`): Python readers for particle and sim-hit ROOT files.
  - JSON (`acts.json`, `acts.examples.json`): geometry and material data, not event data.

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
The [geometry example](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/geometry.py)
shows how to produce such a JSON dump from an existing `TrackingGeometry` (`toFile`/`toJson`),
for example from a source build with DD4hep enabled.
