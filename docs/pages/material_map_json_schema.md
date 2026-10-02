@page material_map_json_schema Material map schema

Version 1 stores **surface material**. @ref Acts::TrackingGeometryMaterialJsonConverter
reads and writes JSON or CBOR, optionally compressed with zstd. `toJson`/`fromJson`
convert documents; `toFile`/`fromFile` handle files. @ref Acts::TrackingGeometryMaterialJsonConverter::Options "Options" controls indentation
and compression level. Applying material to geometry is a separate operation:

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

Geometry IDs use integer components `volume`, `boundary`, `layer`, `passive`,
`sensitive` and `extra`. Omitted components mean **zero, never wildcards**.
`passive` names the shared approach/passive bit field. Packed IDs are not supported.
Stable keys are nonempty and case-sensitive. Duplicate identities are errors;
a proto payload's `material_key`, when present on a keyed assignment, must match.

Absent assignments leave geometry unchanged. Null payloads are allowed only for
unkeyed assignments and are a no-op when applied. Vacuum is an explicit material,
not an absent assignment. Proto materials describe mapping intent; merged markers
record lost material and cannot be applied as physical material.

## Values and units

Format units are **mm**, **radians**, **GeV** and **mol/mm³**, independent of build
units. Relative atomic mass and atomic number are dimensionless.
Material fields preserve radiation and interaction lengths, relative atomic mass,
atomic number, molar density, molar electron density and mean excitation energy.
The last two are stored independently rather than reconstructed from the others.

Vacuum is `{"kind":"vacuum"}`. The string `"infinity"` is allowed for radiation
and interaction lengths; numeric NaN and infinities are outside the format.
Slabs retain thickness, including for vacuum. Surface settings preserve mapping
type (`default`, `pre`, `post`, `sensor`) and split factor in [0,1]. Proto materials
require split factor 1.

## Axes and storage

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

Proto-grid axes may defer ranges, directions and boundary behavior to geometry.
Deferred-variable edges run from 0 to 1 and scale to the resolved range.
Explicit properties must agree with geometry. Resolved material cannot contain
unresolved axes. Proto surface binning may have zero axes for homogeneous mapping.

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
custom kinds also need an extended schema for offline validation.

Legacy conversion follows decode → @ref Acts::TrackingGeometryMaterial "TrackingGeometryMaterial" → new encode.
It cannot recover information already lost by the legacy reader or format,
such as omitted settings, guard cells, store sharing or independent material
properties. Existing legacy converters retain their behavior.
