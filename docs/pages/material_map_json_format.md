@page material_map_json_format Material map JSON format

A material map is the file the material mapping writes and the
@ref Acts::JsonMaterialDecorator "JsonMaterialDecorator" reads back to attach
material to a tracking geometry. This page documents the on-disk layout as
produced by @ref Acts::MaterialMapJsonConverter "MaterialMapJsonConverter" and
@ref Acts::SurfaceMaterialJsonConverter "SurfaceMaterialJsonConverter".

The example below is not typed out by hand. It is generated from the converters
by the `MaterialJsonDocumentation` unit test, which fails if this file stops
matching what the code writes; run that test with `ACTS_UPDATE_DOC_EXAMPLES=1`
to refresh it after a format change. It is deliberately as small as the format
allows, and shows the three surface payloads that differ in shape plus one
volume entry. The remaining payload types differ only in the keys tabulated
further down.

@include examples/material_map_example.json

## Document layout

The document has two top level keys, `Surfaces` and `Volumes`. Each holds a
@ref Acts::GeometryHierarchyMapJsonConverter "geometry hierarchy map" document,
which is a header naming the container -- `acts-geometry-hierarchy-map`, with a
`format-version` and a `value-identifier` -- followed by a flat list of
`entries`.

An entry carries the non-zero levels of its @ref Acts::GeometryIdentifier
(`volume`, `boundary`, `layer`, `approach`, `sensitive`) next to a `value`
object. Levels that are zero are omitted, which is why the homogeneous entry
above shows `volume` and `boundary` but no `layer`. For material maps the
`value` object has a single `material` key. An entry whose `material` is
missing or `null` is skipped when reading, which is how a geometry dump can
list surfaces that carry no material yet.

## Surface material payloads

The value under `material` is a self-describing payload. Its `type` tag selects
the decoder and is authoritative: a missing or unknown tag is an error, the
payload is never guessed from the keys that happen to be present.

| `type`                   | C++ type                                                                                                    | Payload keys                                              |
|--------------------------|-------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------|
| `homogeneous`            | @ref Acts::HomogeneousSurfaceMaterial                                                                        | `data`                                                    |
| `binned`                 | @ref Acts::BinnedSurfaceMaterial                                                                             | `binUtility`, `data`                                      |
| `proto`                  | @ref Acts::ProtoSurfaceMaterial                                                                              | `binUtility`                                              |
| `proto-grid`             | @ref Acts::ProtoGridSurfaceMaterial                                                                          | `axis_specs`                                              |
| `grid`                   | @ref Acts::GridSurfaceMaterial                                                                               | `accessor`                                                |
| `merged-material-marker` | @ref Acts::MergedMaterialMarker                                                                              | none                                                      |

Two more keys are common to all of them:

- `mapMaterial` steers the material mapping. Reading a payload with
  `mapMaterial: false` yields no material at all, so this is also how a surface
  is flagged out of the mapping. Proto material without any binning is written
  with `mapMaterial: false` for that reason.
- `mappingType` is one of `PreMapping`, `Default`, `PostMapping` or `Sensor`
  and tells the mapper where along the propagation the material should be
  assigned. It is absent from the `grid` and `merged-material-marker`
  payloads, which do not participate in the deprecated mapping path.

### Material slabs

Wherever material itself is stored, it is a slab: the opaque
@ref Acts::Material parameter vector plus a thickness. The vector currently
holds radiation length, interaction length, relative atomic mass, nuclear
charge and molar density in ACTS native units, but it is deliberately opaque --
read it through the converter rather than by index, since more parameters may
be appended later. Vacuum is written as a `null` vector, as in the empty second
bin of the binned entry above.

### `homogeneous` and `binned`

`homogeneous` is one slab for the whole surface. The slab sits in a nested
array for historical reasons: it is the degenerate case of the `binned` matrix,
which pairs a @ref Acts::BinUtility with the slabs it addresses. That matrix is
indexed `[bin of the second binning][bin of the first binning]`, so the one
dimensional binning of the example gives a single row of two slabs.

### `proto` and `proto-grid`

Binning instructions for the material mapping that carry no material yet.
`proto` expresses the binning as a @ref Acts::BinUtility, exactly as `binned`
does but without the `data`. `proto-grid` expresses it as an `axis_specs` list
of @ref Acts::AxisSpec, which is the representation the grid based material
uses; each spec has a `type`, a `bins` count, a `range`, a `boundary_type` and
a `direction`. Exactly two specs are required.

### `grid`

The whole grid material family shares one tag, and one C++ class:
@ref Acts::GridSurfaceMaterial. The payload is a grid under `accessor`, made of
a list of `axes` and a `data` list of `[local bins, value]` pairs. The local bin
indices are **one based** and follow the axis order of `axes`.

The grid is always two dimensional, and lookup is local: `loc0` addresses axis
0 and `loc1` axis 1 directly. There is no global (position) lookup and hence no
coordinate-transform description in the payload.

What sits in a bin depends on `accessor.type`, which names the storage backend:

| `accessor.type`    | Bin value       | Extra keys                     |
|--------------------|-----------------|--------------------------------|
| `direct`           | a material slab | none                           |
| `indexed`          | an index        | `storage_vector`               |
| `globally_indexed` | an index        | `storage_vector` **or** `store` |

`indexed`, shown in the example, keeps a slab store next to the grid and the
bins index into it, which pays off as soon as several bins share the same slab.
Note that the store index is unrelated to the bin number: the single bin of the
example holds index 1, the second entry of its `storage_vector`.

`globally_indexed` is the same, except that the store may be shared with other
surfaces. A standalone payload inlines the store as `storage_vector` and stays
self-contained; when a document-wide store table is in use, the entry instead
references it by id under `store`, so the sharing survives the round trip
rather than being flattened into one copy per surface. The two are mutually
exclusive: an entry carrying `store` can only be read with the table that
defines it.

### `merged-material-marker`

A sentinel left behind by @ref Acts::Portal::merge when the material of two
merged portal surfaces had to be dropped. Its payload is just the tag and
`mapMaterial`; it carries no material, only the information that something was
lost here.

## Volume material payloads

Volume material entries are written by a separate, older converter and use
their own tags: `homogeneous` (a single @ref Acts::Material parameter vector
under `data`, as in the example above), `proto` (a `binUtility` only), and
`interpolated2D` / `interpolated3D` (a `binUtility` plus one parameter vector
per bin). The tags overlap with the surface ones by accident; the two lists
live under different top level keys and are never mixed.
