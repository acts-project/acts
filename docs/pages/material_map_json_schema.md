@page material_map_json_schema Material map schema

Version 1 stores **surface material**. @ref Acts::TrackingGeometryMaterialJsonConverter
reads and writes JSON or CBOR, optionally compressed with zstd. `toJson`/`fromJson`
convert documents; `toFile`/`fromFile` handle files. @ref Acts::TrackingGeometryMaterialJsonConverter::Options "Options" controls indentation
compression level and optional material quantization. Applying material to geometry is a separate operation:

@snippet{trimleft} examples/material_map_json.cpp Read and write material map

Python exposes the converter as `acts.json.TrackingGeometryMaterialJsonConverter`.
Volume material remains available through the legacy format; version 1 rejects it.

The schema is in `Plugins/Json/schema/material-map-v1.schema.json` and installs
under `share/Acts/schema`. It defines fields, variants and numeric bounds.
Examples in `docs/examples/material-map-v1/` cover homogeneous material
(`minimal.json`), mapped surfaces (`surfaces.json`) and mapping templates
(`templates.json`). Their relative `$schema` references support editor validation.

@include examples/material-map-v1/minimal.json

## Document and assignments

`header` contains `format: "acts-material-map"`, `version: 1`, and an optional
human-readable `description`. Unknown versions are rejected. Read or change the
description through @ref Acts::TrackingGeometryMaterial::description "description()" and @ref Acts::TrackingGeometryMaterial::setDescription "setDescription()".
The optional root `$schema` identifies the schema; it is not fetched by the reader.

`surfaces` contains material assignments; `slab_stores` holds shared grid data.
An empty assignment list is valid.

| Target kind | Matching |
| --- | --- |
| `geometry-id` | Exact `geometry_id`, without hierarchy wildcards |
| `stable-key` | Exact `key`; `recorded_geometry_id` is diagnostic only |

Geometry IDs use integer components `volume`, `portal`, `layer`, `passive`,
`sensitive` and `extra`. Omitted components mean **zero, never wildcards**.
`portal` names the boundary bit field; `passive` names the shared approach/passive
bit field. Packed IDs are not supported.
Stable keys are nonempty and case-sensitive. Duplicate identities are errors;
a proto payload's `material_key`, when present on a keyed assignment, must match.

Absent assignments leave geometry unchanged. Null payloads are allowed only for
unkeyed assignments and are a no-op when applied. Vacuum is an explicit material,
not an absent assignment. Proto materials describe mapping intent; merged markers
record lost material and cannot be applied as physical material.

## Values and units

The required `header.units` declares base units by name. Writers use:

```json
"units": {"length": "mm", "angle": "rad", "energy": "GeV", "material_amount": "mol"}
```

| Quantity | Accepted units |
| --- | --- |
| `length` | `nm`, `um`, `mm`, `cm`, `m` |
| `angle` | `rad`, `mrad`, `deg` |
| `energy` | `eV`, `keV`, `MeV`, `GeV`, `TeV` |
| `material_amount` | `mol`, `mmol`, `kmol` |

Readers convert to Acts units and reject unknown or dimensionally wrong names.
Both molar densities use **material_amount / length³**; there is no compound-unit
parser. Changing a length unit therefore also requires rescaling density values.
Lengths cover slab thickness, material lengths, translations and length-valued
axes. `phi`/`theta` use the angle unit; `eta` and normalized deferred edges are
dimensionless. Axis ranges without a direction use the length unit; `rphi` uses
length with azimuth in radians. Relative atomic mass and atomic number are
dimensionless.
Material fields preserve radiation and interaction lengths, relative atomic mass,
atomic number, molar density, molar electron density and mean excitation energy.
The last two are stored independently rather than reconstructed from the others.

Vacuum is `{"kind":"vacuum"}`. The string `"infinity"` is allowed for radiation
and interaction lengths; numeric NaN and infinities are outside the format.
Slabs retain thickness, including for vacuum. Surface settings preserve mapping
type (`default`, `pre`, `post`, `sensor`) and split factor in [0,1]. Proto materials
require split factor 1.

## Optional quantization

@ref Acts::TrackingGeometryMaterialJsonConverter::Options::materialFractionBits
controls the retained float32 fraction bits (0–23); the default 23 leaves values
unchanged. Lower values improve compression by rounding built-in material
properties and slab thickness to fewer bits. Rounding uses nearest, ties to even,
with relative error at most `2^(-bits-1)` for normal values. Zeros, subnormals,
infinities and values that would overflow are preserved.

Geometry, axes, settings and custom payloads are unchanged. Quantization applies
to canonical output units and does not change the schema or reader. It is
available through `toJson` and `toFile`, including the Python options.

## Axes and storage

@ref Acts::BinUtility "BinUtility" is being phased out. Its binned/proto and
subdivided-axis encodings are retained for compatibility. Prefer
@ref Acts::GridSurfaceMaterial and @ref Acts::ProtoSurfaceMaterial for new code.

Arrays are dense, flat and zero-based, with **axis 0 varying fastest**:
`offset = i0 + size0 * i1`. This differs from native @ref Acts::MultiAxis "MultiAxis" storage order;
the converter handles the permutation.

- Binned surfaces store regular cells only, with one or two axes.
- Grid surfaces have two resolved axes and include underflow and overflow cells
  on both axes, including bound/closed axes. Each extent is `bins + 2`.
- `open` uses guard cells, `bound` clamps, and `closed` wraps. Intervals include
  the lower edge and exclude the upper edge. @ref Acts::BinUtility "BinUtility" `open` maps to `bound`.
- Equidistant axes have a range and bin count; variable axes have increasing
  edges. Grid lookup uses surface-local coordinates in axis order.
- Binned/proto material may carry a rigid local-to-global transform: a row-major
  3×3 rotation and a translation vector. An omitted transform is identity.
- Subdivided @ref Acts::BinUtility "BinUtility" axes preserve their base and subdivision: `replace`
  refines one matching interval; `repeat` tiles the subdivision in every base bin.
  Legacy `theta`/`mag` @ref Acts::BinUtility "BinUtility" directions and custom coordinate callbacks
  are not supported.

Writers emit `proto-grid` for @ref Acts::ProtoSurfaceMaterial, whose two axes
may defer ranges, directions and boundary behavior to geometry.
Deferred-variable edges run from 0 to 1 and scale to the resolved range.
Explicit properties must agree with geometry. Resolved material cannot contain
unresolved axes. Proto surface binning uses two one-bin axes for homogeneous
mapping.

| Grid storage | Slab data |
| --- | --- |
| `direct` | Slab values in each cell |
| `indexed` | Cell indices into a local `slabs` list |
| `globally-indexed` | Cell indices into a named entry in `slab_stores` |

Indices must be in range. References to the same store retain shared allocation;
separate stores remain separate even if their values are equal.

## Validation and compatibility

Run `python CI/check_material_schema.py` with `jsonschema>=4.18` to validate the
schema and examples, or use the pre-commit hook. Schema tests run with
`python -m pytest CI/check_material_schema.py`.

The schema checks structure and physical-value bounds. The C++ reader checks
format/version, sizes, identities and references needed to construct material;
it does not perform full schema validation. Geometry application checks matching
and resolves deferred axes. Schema validity alone does not guarantee that a map
can be applied to a particular geometry.

The schema rejects unknown fields. The reader ignores unused fields except
`volumes`, which is rejected. @ref Acts::TrackingGeometryMaterialJsonConverter::Config "Config" supports custom material dispatchers;
custom kinds also need an extended schema for offline validation. Encoder and
decoder contexts expose unit factors for custom payloads; encoding divides by
the factor and decoding multiplies by it.

Legacy conversion follows decode → @ref Acts::TrackingGeometryMaterial "TrackingGeometryMaterial" → new encode.
It cannot recover information already lost by the legacy reader or format,
such as omitted settings, guard cells, store sharing or independent material
properties. Existing legacy converters retain their behavior.
