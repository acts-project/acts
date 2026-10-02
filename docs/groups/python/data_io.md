@defgroup python_data_io Data Reading and Writing
@ingroup python_bindings
@brief Readers and writers available in the PyPI wheel and optional source builds.

# Event input and output

Add readers with `s.addReader(...)`, converters with `s.addAlgorithm(...)`, and writers with
`s.addWriter(...)`. The available formats depend on the build (see @ref python_bindings):

| Format | Python module | Availability and use |
| :--- | :--- | :--- |
| CSV | `acts.examples` | PyPI or source; readers and writers for common ACTS event collections. |
| Parquet | `acts.examples.arrow` | PyPI or source; reads and writes Arrow tables. Converters connect tables to ACTS event collections. |
| HepMC3 | `acts.examples.hepmc3` | PyPI or source; ASCII event files without ROOT. ROOT-format HepMC3 files require ROOT support. |
| ROOT via Uproot | `acts.examples.uproot` | PyPI or source; Python readers for ACTS particle and sim-hit ROOT files. Install `uproot` and `numpy` separately. |
| ROOT | `acts.examples.root` | ROOT-enabled source build; native readers and writers for particles, sim hits, tracks, vertices, material, and performance output. |
| EDM4hep/podio | `acts.examples.edm4hep` | EDM4hep-enabled source build; `PodioReader` and `PodioWriter`. |

## Parquet and Arrow tables

There is one `ParquetReader` for the configured Parquet collections, rather than a separate
reader for each ACTS data type. It places one Arrow table per collection on the event whiteboard.
Configure a directory and an expected schema for each collection:

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

The example makes an Arrow table available as `particles_arrow`; it does not create an ACTS
`SimParticleContainer`. Add an input converter for the dataset's schema when ACTS algorithms need
event collections. The current bindings include `ColliderMLRelease1InputConverter` for the
ColliderML Release 1 particle, hit, and optional track tables; other Parquet schemas need a
matching converter.

For output, `ArrowParticleOutputConverter`, `ArrowSimHitOutputConverter`, and
`ArrowTrackOutputConverter` turn ACTS collections into Arrow tables. One `ParquetWriter` can then
write those tables to separate collection directories. See the
[Parquet round-trip test](https://github.com/acts-project/acts/blob/main/Python/Examples/tests/test_arrow.py)
for a particle example and the
[full-chain example](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/full_chain_odd.py)
for particle, hit, and track output.

## HepMC3 event files

`HepMC3Reader` reads ASCII `.hepmc`/`.hepmc3` files into HepMC3 events. Use
`HepMC3InputConverter` to turn those events into ACTS particles and vertices. In the other
direction, `HepMC3OutputConverter` creates HepMC3 events for `HepMC3Writer`. These ASCII paths
work without ROOT. Reading or writing HepMC3 ROOT files requires a ROOT-enabled build with
HepMC3 ROOT I/O support.

## Reading ACTS ROOT files with Uproot

`acts.examples.uproot` provides `UprootParticleReader` and `UprootSimHitReader`: pure-Python
`IReader`s that read the exact file format written by `RootParticleWriter`/`RootSimHitWriter`,
using [uproot](https://uproot.readthedocs.io/) without a ROOT-enabled ACTS build. They require
`uproot` and `numpy` to be installed separately. For their `IReader` implementation, see
@ref python_custom_algorithms.

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
