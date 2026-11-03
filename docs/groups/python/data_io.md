@defgroup python_data_io Data Reading and Writing
@ingroup python_bindings
@brief Readers and writers available in the PyPI wheel and optional source builds.

# Event data input and output

Add readers with `s.addReader(...)`, converters with `s.addAlgorithm(...)`, and writers with
`s.addWriter(...)`. The available formats depend on the build (see @ref python_bindings):

| Format | Python module | Available in PyPI | Use |
| :--- | :--- | :---: | :--- |
| CSV | `acts.examples` | ✓ | Readers and writers for common ACTS event collections. |
| Parquet | `acts.examples.arrow` | ✓ | Reads and writes Arrow tables; converters connect them to ACTS event collections. |
| HepMC3 (ASCII) | `acts.examples.hepmc3` | ✓ | Reads and writes HepMC3 event files without ROOT. |
| ROOT (native) | `acts.examples.root` | ✗ | Readers and writers for particles, sim hits, tracks, vertices, material, and performance output. |
| ROOT via Uproot | `acts.examples.uproot` | ✓ | Python readers for ACTS particle and sim-hit ROOT files; install `uproot` and `numpy` separately. |
| EDM4hep/podio | `acts.examples.edm4hep` | ✗ | `PodioReader` and `PodioWriter`. |

✓ is included in the PyPI wheel; ✗ requires a source build with the corresponding component
enabled. HepMC3 ROOT-format files also require ROOT support.

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

# Loading detector geometry from JSON

To use a `TrackingGeometry` built with detector tooling such as DD4hep or TGeo in a PyPI-only
script, export it to JSON in a source build and load the JSON file:

```python
import acts.json

gctx = acts.GeometryContext()
converter = acts.json.TrackingGeometryJsonConverter()
trackingGeometry = converter.fromFile(gctx, "geometry.json")
```

Material maps can be attached the same way with `acts.json.JsonMaterialDecorator`.
The [geometry example](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/geometry.py)
shows how to produce such a JSON dump from an existing `TrackingGeometry` (`toFile`/`toJson`),
using a source build with the relevant detector tooling enabled.
