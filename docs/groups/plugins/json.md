@defgroup json_plugin JSON Plugin
@ingroup plugins
@brief JSON (de)serialization of geometry, material and configuration.

The file format written by the material converters is documented in
@ref material_map_json_format.

The @ref material_map_json_schema describes the versioned material
format handled by @ref Acts::TrackingGeometryMaterialJsonConverter. Legacy material
converters are deprecated and retain the previous format.

## Tracking geometry with bundled material

`TrackingGeometryJsonConverter` writes geometry format version 2. The geometry
header contains `format-version: 2`; version 1 geometry with inline surface
material remains readable. Version 2 never writes inline material payloads.

By default, the top-level `material` member contains a complete versioned material
document, produced by `TrackingGeometryMaterialJsonConverter`. Its format version
is independent of the geometry version. Material stores and application-defined
material kinds use the same codec as standalone material files. Pass the material
converter's full `Config` as `Config::materialConfig`, and its encoding `Options`
as `Options::materialOptions`. The outer geometry options control file indentation
and zstd compression; nested material options control material precision.

A surface with a keyed proto-material carries `material_key` in its geometry
record. This identity remains present when `Options::writeMaterial` is false.
Reading restores a keyed proto-material placeholder and applies the bundled
material with the ordinary material assignment rules: stable keys for keyed
surfaces, geometry IDs otherwise. Missing keys and incompatible assignments fail
as they do for standalone material application; unused assignments are permitted.
Version 2 does not regenerate surface geometry IDs during import.

Applying mapped material replaces proto-material placeholders. Extracting material
from an already decorated geometry therefore produces geometry-ID assignments;
the bundle is a snapshot and does not recover the original mapper output's stable
keys. Keep that original keyed material document for use with rebuilt geometries.

Set `Options::materialOptions.materialFractionBits` below 23 to round material
properties and thickness; geometry coordinates are unaffected. Unsupported material,
including volume material, is rejected by the material codec. A geometry-only
export omits the material document entirely.

Use `toFile` and `fromFile` for file I/O. A `.json.zst` filename selects JSON with
zstd compression; input encoding and compression are detected from the content.
