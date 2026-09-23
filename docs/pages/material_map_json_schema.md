@page material_map_json_schema Material map schema proposal

**Experimental version 1.** `Acts::TrackingGeometryMaterialJsonConverter` reads
and writes this format. Its `Config` supplies extensible surface dispatchers;
`toJson`/`fromJson` handle documents and `toFile`/`fromFile` handle files. `Options`
controls indentation, compression level and optional material quantization. Applying material remains a separate `TrackingGeometryMaterial::apply`
operation. Existing material writers, decorators and converters are deprecated
and retain their legacy-format behavior; the new codec does not call them.
The schema remains a draft pending review of the supported-state boundaries below.
The hand-authored examples are exercised by C++ codec tests as well as offline QA.

The schema is shipped in `Plugins/Json/schema/material-map-v1.schema.json`.
Examples in `docs/examples/material-map-v1/` use a relative `$schema` reference for
editor completion and offline validation. The schema's URN is an identifier, not
a download endpoint. No network resolution or C++ schema-validator dependency is
needed. Run `python CI/check_material_schema.py` with Python `jsonschema>=4.18`
installed; the pre-commit hook supplies that dependency in its own environment.
The same file contains pytest self-tests, run separately by the CI tooling
self-test job with `python -m pytest CI/check_material_schema.py`.
The new API is declared in
`ActsPlugins/Json/TrackingGeometryMaterialJsonConverter.hpp`. For example:

```cpp
Acts::TrackingGeometryMaterialJsonConverter converter;
auto material = converter.fromFile("material.cbor.zst");
material.apply(geometry);
converter.toFile(material, "material.json");
```

The JSON plugin also installs `ActsMaterialMapMigrate` in `bin` to migrate
files written by the deprecated material converter:

```sh
ActsMaterialMapMigrate old-material.json material.json
ActsMaterialMapMigrate old-material.json.zst material.cbor.zst \
    --material-fraction-bits 16 --compression-level 19
```

The input encoding is detected from its contents. The output extension selects
JSON or CBOR, with optional zstd compression. Defaults preserve full float32
precision, use four-space indentation and zstd level 9. `--indentation` changes
text indentation; `--help` lists all options. Volume material causes migration
to fail because version 1 only supports surfaces. The tool preserves material
assignments and stable keys supported by the legacy reader; unrelated geometry
annotations in decorated legacy files are not part of the new material format.
Migration uses the legacy reader's semantics, including its default split factors
(the legacy format does not store them) and normalization of single-bin axes;
it cannot recover settings already lost by the legacy format.

Material quantization is opt-in and affects only built-in slab thickness and
composition fields. `Options::materialFractionBits` defaults to 23 (full float32
precision); values from 0 to 22 round to fewer binary fraction bits. For example:

```cpp
Acts::TrackingGeometryMaterialJsonConverter::Options options;
options.materialFractionBits = 16;
converter.toFile(material, "material.cbor.zst", options);
auto document = converter.toJson(material, options);
```

Rounding uses nearest with ties to even. For `b` retained fraction bits, each
normal value changes by at most `2^(-b-1)` relative to its original value (about
0.000763% at 16 bits). Signed zero, subnormal values, infinity sentinels and values
that would round to infinity are preserved. This bound applies to individual
stored values, not to derived physics quantities or tracking results. Repeated
serialization at the same precision is idempotent. The input material stays
unchanged, and the reader needs no quantization setting or schema change.
Geometry coordinates, settings, identifiers and custom surface payloads retain
their existing precision. Shared slab stores use the same quantization.
Binary quantization can improve compressibility but does not necessarily produce
short decimal text; CBOR continues to store the resulting values as float32.

Plugin installation also places the schema under `share/Acts/schema` (subject
to `CMAKE_INSTALL_DATADIR`). The schema remains explicitly marked as a draft.

| Example | Coverage |
| --- | --- |
| `minimal.json` | One homogeneous surface |
| `surfaces.json` | Binned/direct/indexed/shared grids, both identities, description, vacuum and null |
| `templates.json` | Proto surfaces, deferred proto-grid, nested binning, merge origins |

@include examples/material-map-v1/minimal.json

## Envelope and identity

`format: "acts-material-map"` and `version: 1` govern the entire document,
including axes, assignments, materials and stores. There are no nested versions.
`$schema` is optional tooling information, never permission to bypass version
checks. JSON text and CBOR encode the same data model; zstd only compresses it.
CBOR must use the same string keys, IDs and finite numeric values as JSON.
Unknown versions are rejected. The schema rejects unknown fields; the runtime
reader ignores unused fields except the explicitly unsupported `volumes` field.
The shipped schema enumerates the built-in variants and rejects custom kinds.
An application may register additional material encoders/decoders with the type
dispatcher: its reader/writer can support those kinds even though the standard
schema rejects them. Offline validation of such documents requires an extended
schema supplied by the application. An unregistered kind must fail clearly at
runtime; optional schema validation does not disable dispatcher extensibility.

Version 1 supports **surface material only**. `surfaces` is the assignment list;
an empty list is valid. There is no `volumes` field: readers reject it, and writers
reject any nonempty `TrackingGeometryMaterial::volumeMaterials`, including null
assignments. Volume material remains available through the legacy format; a
future document version can add support without changing version 1's meaning.
Payloads live on assignments; `slab_stores` shares slab data between surface
grids. Stable keys and store names are nonempty, case-sensitive strings. No
normalization, path interpretation or external reference fetching is implied.

| Assignment | Identity | Other information |
| --- | --- | --- |
| Surface `geometry-id` | Exact `geometry_id` | No hierarchy wildcard lookup |
| Surface `stable-key` | Exact `key` | Required `recorded_geometry_id` is diagnostic only |

IDs are readable objects with integer components: `volume` (0–255), `boundary`
(0–255), `layer` (0–4095), `approach` (0–255), `sensitive` (0–1048575), and `extra`
(0–255). These cover the complete current uint64 identifier without large JSON
numbers. Omitted components mean **zero, never wildcards**; `{}` is the zero ID,
including pre-closure provenance. Writers omit zero components; readers also
accept explicit zeros. `approach` is the canonical spelling of the shared
approach/passive field; there is no separate `passive` property. The same object
form is used for assignments, recorded geometry IDs and merge origins.
Packed integer and decimal-string forms are not part of this draft.
Duplicate identities within each assignment namespace are errors. Different keys
may record the same ID, and keyed/unkeyed entries may coexist. Never fall back from
a stable key to the recorded ID. An optional proto `material_key` belongs to the
payload (including unkeyed designators); if its assignment is keyed, they must
agree. Merge origins are provenance, not executable assignments or references.

An absent assignment leaves geometry unchanged. A null unkeyed surface payload
preserves a null container entry and is a no-op under current `apply`. Keyed payloads
must be non-null. Neither null nor absence means vacuum: a homogeneous vacuum
material is an actual assignment. A slab containing vacuum with thickness zero is
`MaterialSlab::Nothing()`; vacuum with positive thickness retains that thickness.
An ordinary material with zero slab thickness also retains its composition.
Proto materials describe mapping intent and are never discarded based on a
`mapMaterial` flag. A merged marker records lost material and is not usable
physical material; geometry application retains the existing marker rejection.

## Physical values and units

All lengths and translations are **mm**, angles **radians**, energies **GeV**,
and molar densities **mol/mm³**. Relative atomic mass and atomic number are
dimensionless. These are fixed format units, independent of a build's units;
C++ conversion uses `UnitConstants`. `theta` is a polar angle and
`eta` is pseudorapidity; `rphi` is transverse radius times azimuth, and `mag` is
3D radius. Axis direction determines its unit. No configurable unit header.

A material has seven explicit physical fields. Read/write via `X0`, `L0`, `Ar`,
`Z`, `molarDensity`, `molarElectronDensity`, `meanExcitationEnergy` and the full
`fromMolarDensity` constructor, with unit conversion. This is a deliberate draft
exception to `Material.hpp`'s advice to serialize its opaque parameter vector:
the current five-component vector loses independently supplied electron density
and excitation energy. Do not assume a vector's order or length in the new codec.
No mass-density conversion or automatic recomputation of the two extra fields.
Float conversion requires range checks; bitwise equality after decimal/unit
conversion is not promised. Physical-value constraints belong to the schema;
the codec does not duplicate them or revalidate encoded output.

`"infinity"` is allowed only for radiation and interaction lengths. Vacuum is
`{"kind":"vacuum"}`, not a collection of nonfinite numbers. JSON/CBOR NaN and
numeric infinities are forbidden everywhere. Slab thickness is nonnegative.
Surface settings explicitly preserve mapping type and pre/post split factor:
`default`, `pre`, `post`, `sensor` correspond to the four `MappingType` values;
1 favors opposite-pre and 0 along-pre. The draft restricts the factor to [0,1].

## Binning, ordering and coordinates

All data arrays are dense and flat. **Axis 0 varies fastest**:
`offset = i0 + size0 * (i1 + size1 * i2)`. All serialized positions and store
indices are zero-based. There are no implicit holes, default index zero, skipped
null entries, duplicate cell coordinates or inferred dimensions.
This differs from the current MultiAxis storage, where the trailing axis varies
fastest: a codec must explicitly permute grid storage, not copy its flat
vector. The common external order follows the existing binned surface matrix.

* Binned surface arrays contain regular cells only, with product-of-bin-counts
  entries and one or two axes. Legacy matrices `[i1][i0]` flatten in this order.
* Grid surfaces have exactly two resolved axes. Arrays include a full guard
  shell, even for bound/closed axes: each extent is `bins + 2`. Position 0 is
  underflow, 1 through N are regular cells, N+1 is overflow. This explicitly
  preserves storage accepted by the full `GridSurfaceMaterial` constructor,
  beyond the regular-only convenience constructors and legacy JSON output.
* `open` routes out-of-range coordinates to guard cells, `bound` clamps to the
  first/last regular cell, and `closed` wraps periodically. Regular intervals are
  lower-inclusive/upper-exclusive; the top edge therefore overflows, clamps or
  wraps respectively. Guard cells on bound/closed axes are retained even if
  normal lookup cannot reach them. BinUtility `open` maps to **bound**, not open.
* Equidistant axes require increasing range endpoints and positive counts;
  variable edges strictly increase. Grid directions are optional annotations:
  lookup consumes surface-local coordinates directly in axis order, without a
  transform or another projection. Geometry must supply compatible coordinates.
  BinUtility axes require directions and only bound/closed behavior, since their
  storage has no guard cells.
* Binned/proto `binning` retains an optional rigid **local-to-global** transform.
  For global lookup, invert it first, then extract directions in the stated order.
  Surface-local lookup already consumes the appropriate axis coordinates.
  Rotation is a 3×3 row-major proper orthonormal matrix, translation a 3-vector;
  omitted transform is identity. Arbitrary executable coordinate callbacks
  cannot be serialized and must cause the writer to reject the object.
* `subdivided` axes preserve BinUtility's nested mapping intent: `replace`
  substitutes the subdivision into the single base interval matching its range;
  `repeat` tiles its offsets into every equal-width base interval. Direction and
  boundary behavior come from the base; nested axes must agree. Repeat requires
  a subdivision range equal to the first base interval. Counts are computed after
  recursive refinement. Resolved grids use flattened edges instead.

Proto surface settings require split factor 1, matching the Core proxy types.
Proto surface binning permits zero axes (homogeneous mapping intent).
Proto-grid has exactly two axis specs: equidistant range, boundary and direction
may be omitted independently; variable edges fix a range but may defer boundary
and direction. Deferred-variable edges strictly increase from exactly 0 to 1;
resolution maps them to `min + edge * (max - min)`. Missing properties come from
the geometry; supplied properties must agree with it, not override it. Deferred
specs are not allowed in resolved payloads. `templates.json` includes homogeneous,
nested and deferred templates and a merged marker.

Direct grid storage holds slabs; indexed storage owns a local `slabs` list;
globally-indexed storage references one document `slab_stores` entry. All indices
must be in range. Repeated references to one name retain shared allocation in a
decoder; distinct names do not imply sharing even when contents compare
equal. There is no global-store inline fallback. Standalone export must carry its
store table. Unused stores are legal. Mutation/scaling semantics remain those of
the material container, not a promise of independent per-assignment copies.

## Description and validation boundaries

An optional top-level `description` string provides human-readable context. It
has no effect on material application. Read it with `material.description()`;
set it with `material.setDescription("Run A")`. Passing `std::nullopt` removes
the field; an empty string is preserved. No general metadata, provenance log,
or extension container is supported. Format/version/schema reference remain
serializer concerns.

| Layer | Responsibility |
| --- | --- |
| JSON Schema | Required fields, types, disjoint variants, fixed dimensions, basic numeric bounds, no unexpected structural fields |
| C++ document reader | Format/version, duplicate JSON keys, integer and allocation bounds, array sizes, unique identities, key consistency, store references and indices, and unsupported representations. Required fields and types use JSON library exceptions; axis constructors validate axis parameters. |
| Geometry application/resolution | Existing key matching, uniqueness of participating surfaces, Gen1/key incompatibility, marker rejection, coordinate/bounds agreement and deferred-axis resolution; description is ignored |

Schema validity alone never proves a map is physically meaningful or applicable
to a detector. `CI/check_material_schema.py` validates examples structurally and
checks representative document invariants as fixture QA; it is not a public
material reader or complete semantic validator. Use the offline schema check
for the full format contract. The runtime codec performs conversion and basic
construction checks, not exhaustive schema or physical-value validation. File decoding rejects duplicate keys before constructing
the JSON/CBOR DOM; `fromJson` cannot recover duplicates already overwritten by a
caller's parser. File documents are limited to 256 nesting levels and recursive binning to 32.

## Legacy migration and remaining review decisions

The intended flow is legacy decode → `TrackingGeometryMaterial` → new encode.
Before that flow can be described as lossless, migration must inspect and report
legacy losses: `mapMaterial:false` drops proto/other payloads; `subdata` is written
but not reconstructed; split factors and grid mapping settings are omitted;
regular-only grids lose boundary cells; missing/invalid grid entries get defaulted;
standalone global stores lose sharing;
opaque five-number materials already lost electron/excitation overrides. Once
lost by decoding, new serialization cannot recover them. Migration must either
reject affected input or use an explicitly reported preservation-aware legacy
path. Never infer lost values and call the result lossless.

The schema intentionally makes stricter choices (physical values, split
range, complete cells) than all states C++ can construct.
Review before freezing: legacy `theta`/`mag` BinUtility projection quirks. These two BinUtility directions
are rejected rather than being reinterpreted as their mathematical names.
Nested binning retains the original base/subdivision construction rather than
flattening its intent. Existing legacy paths remain available with deprecation
warnings; no ROOT cleanup or geometry-format migration is included.
