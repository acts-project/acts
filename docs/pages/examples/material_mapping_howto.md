@page material_mapping_howto How to map material for your detector

# How to map material for your detector

This page walks through producing a material map for a detector of your own. It
assumes you already have the detector described both as an ACTS
@ref Acts::TrackingGeometry and as a detailed simulation geometry that Geant4
can navigate (typically DD4hep, which gives you both from one source).

For the conceptual picture of what the mapping does, see @ref material_mapping.
For the reference description of the Examples chain and the tests that guard it,
see @ref material_mapping_workflow.

The chain is:

```
geometry  ──▶  choose surfaces  ──▶  record (Geant4)  ──▶  map  ──▶  use  ──▶  validate
```

Only the "choose surfaces" step is detector-specific work. The rest is running
three scripts.

## Step 1: choose which surfaces carry material

Mapping only writes material onto surfaces you have explicitly nominated. A
fresh detector nominates nothing, so this is where you start.

### 1a. Dump the geometry to JSON

```console
python Examples/Scripts/Python/geometry.py
```

This writes `geometry-map.json` in the current directory (among other outputs).
The relevant part is in `runGeometry()` in
`Examples/Scripts/Python/geometry.py`, which configures a
`MaterialMapJsonConverter` and a `JsonMaterialWriter`.

> [!important]
> The converter must be configured with `processNonMaterial=True`. Surfaces that
> do not already carry material are otherwise skipped entirely
> (`Plugins/Json/src/MaterialMapJsonConverter.cpp:383`), so they never appear in
> the dump and you have nothing to switch on. `geometry.py` already sets this;
> if you write your own dumping script, do not omit it.

Point the script at your own detector by replacing the `getOpenDataDetector()`
call at the bottom of the file.

### 1b. Reduce it to an editable config

`geometry-map.json` has one entry per surface, which for a real detector is
thousands of entries. `writeMapConfig.py` collapses it into one representative
entry per surface *category* per volume:

```console
python Examples/Scripts/MaterialMapping/writeMapConfig.py geometry-map.json config-map.json
```

### 1c. Edit the config

In `config-map.json`, for every surface category you want to map:

- set `"mapMaterial": true`
- set the `"bins"` entries under `"binUtility" -> "binningdata"` to the
  granularity you want

A homogeneous surface is `1 x 1`. Start coarse: bins with no hits produce no
material, and a fine binning with too few recorded tracks gives you a map full
of holes.

### 1d. Write the choices back

```console
python Examples/Scripts/MaterialMapping/configureMap.py geometry-map.json config-map.json
```

> [!warning]
> This rewrites `geometry-map.json` **in place**. Keep a copy of the pristine
> dump, or be ready to regenerate it with `geometry.py`.
>
> It also only copies the bin *counts* back, not the bin ranges or axis types —
> those come from the geometry dump. Editing `"min"`, `"max"` or `"type"` in
> `config-map.json` has no effect.

`geometry-map.json` is now your **mapping configuration**: the same file you
pass as `--matconfig` in step 3.

## Step 2: record the material with Geant4

```console
python Examples/Scripts/Python/material_recording.py -n 1000 -t 1000 -o mydet_geant4
```

Shoots geantinos through the *detailed* geometry and records what they traverse,
into `mydet_geant4.root`. `-n` is events, `-t` is tracks per event, so the above
is a million tracks; `--eta-range` and `--phi-range` restrict the solid angle.

This step is independent of your surface choices — you only need to redo it if
the detailed geometry changes, not when you retune binning. It is by far the
slowest step, so record generously once and reuse the file.

## Step 3: run the mapping

```console
python Examples/Scripts/Python/material_mapping.py \
    -n 1000000 -i mydet_geant4.root --matconfig geometry-map.json -o mydet_material
```

`--matconfig` loads your configuration through
`acts.IMaterialDecorator.fromFile()`, which puts an
@ref Acts::ProtoSurfaceMaterial (a binning specification with no content yet) on
every surface you enabled. The script then reads those surfaces back out with
`trackingGeometry.extractMaterialSurfaces()` and hands the resulting list to
@ref Acts::IntersectionMaterialAssigner and
@ref Acts::BinnedSurfaceMaterialAccumulator.

Outputs:

| File | Contents |
|---|---|
| `mydet_material_map.json` | the material map, human-readable |
| `mydet_material_map.root` | the same map, for production use |
| `mydet_material_mapped.root` | recorded interactions that found a surface |
| `mydet_material_unmapped.root` | recorded interactions that did not |

`_unmapped.root` is the one to look at when something is wrong. A large unmapped
fraction means the material had nowhere to go — see
@ref material_mapping_howto_troubleshooting "Troubleshooting".

## Step 4: use the map

```python
decorator = acts.IMaterialDecorator.fromFile("mydet_material_map.root")
detector = getMyDetector(materialDecorator=decorator)
```

`.json`, `.cbor` and `.root` are all accepted. This is exactly how the ODD picks
up `data/odd-material-maps.root` by default.

## Step 5: validate

```console
python Examples/Scripts/Python/material_validation.py \
    -n 1000 -t 1000 -m mydet_material_map.root -o mydet_validated -p
```

This re-records material, now from your mapped map instead of from Geant4, so
you can compare the two. `-p` additionally runs a real propagator with a
navigator and writes `mydet_validated_propagated.root`.

Comparing the default and `-p` outputs is worth doing: the default collection is
navigation-independent, so a difference between them is a *navigation* problem
(material the navigator does not reach), not a mapping problem.

To compare against the Geant4 input:

```console
python Examples/Scripts/MaterialMapping/material_comparison.py
root -l Examples/Scripts/MaterialMapping/Mat_map.C
```

`Examples/Scripts/MaterialMapping/material_mapping_check.py -i mydet_material_mapped.root`
plots how far each interaction was moved to reach its assigned surface, which is
the quickest way to spot material being attached to the wrong thing.

## Tuning

Iterating on binning does **not** require re-recording. Edit `config-map.json`,
re-run `configureMap.py` and step 3 against the same recorded file.

Rules of thumb:

- Bins that receive no tracks stay empty. If your map has holes, either coarsen
  the binning or record more tracks.
- Boundary and approach surfaces generally need less granularity than sensitive
  layers.
- Check `_unmapped.root` after every change, not just at the end.

@anchor material_mapping_howto_limits
## Limitations you should know about

These are properties of the current navigation-less mapper, not of your setup.
Both are cases where the legacy propagation-based mapper did more.

**Volume material is not produced.** @ref Acts::MaterialMapper retrieves volume
assignments from the assignment finder and then discards them; only surface
material is accumulated, and `finalizeMaps()` returns an empty volume map
(`Core/src/Material/MaterialMapper.cpp`). If you need volume material, the
deprecated @ref Acts::VolumeMaterialMapper is currently the only path.

**`mappingType` is not honoured.** The `"mappingType"` key round-trips through
the JSON and is stored on the material, but the current assignment does plain
nearest-intersection matching — its own comment says *"no pre/post matching"*
(`Core/src/Material/MaterialInteractionAssignment.cpp`). `PreMapping`,
`PostMapping` and `Sensor` are only interpreted by the deprecated
@ref Acts::SurfaceMaterialMapper. To steer assignment today, use the
`globalVetos`, `localVetos` and `reAssignments` hooks in
@ref Acts::MaterialInteractionAssignment::Options.

@anchor material_mapping_howto_troubleshooting
## Troubleshooting

**Everything lands in `_unmapped.root`.** Nothing was nominated. Confirm that
`geometry-map.json` actually contains `"mapMaterial": true` somewhere — if you
forgot `processNonMaterial=True` in step 1a, or forgot to run `configureMap.py`
in step 1d, the file will be syntactically fine and semantically empty.

**The map is full of empty bins.** Too fine a binning for the number of recorded
tracks. Coarsen it, or record more.

**Material appears on the wrong surface.** Expected when candidate surfaces are
close together: assignment picks the nearest intersection with no notion of
"before" or "after". Use `material_mapping_check.py` to see the assignment
distances, and the veto/re-assignment hooks to correct specific surfaces.

**Validation disagrees with Geant4, but only with `-p`.** That is navigation,
not mapping. The mapped material is fine; the navigator is not finding all of
it.

## Worked example

`Python/Examples/tests/test_material_mapping.py` runs the whole chain against
the ODD and is the most reliable executable reference. The ODD itself ships a
finished map at `data/odd-material-maps.root`.
