"""Map and validate binned/grid material with deterministic material tracks."""

from array import array
import json
import math

import pytest

import acts
from acts.json import TrackingGeometryMaterialJsonConverter


def material_surface(tmp_path):
    surface = acts.Surface.createPlane(
        acts.Transform3.Identity(), acts.RectangleBounds(10.0, 10.0)
    )
    surface.assignGeometryId(acts.GeometryIdentifier(sensitive=1))
    document = {
        "Surfaces": {
            "acts-geometry-hierarchy-map": {
                "format-version": 0,
                "value-identifier": "Material Surface Map",
            },
            "entries": [
                {
                    "sensitive": 1,
                    "value": {
                        "material": {
                            "type": "proto-grid",
                            "mapMaterial": True,
                            "mappingType": "Default",
                            "axis_specs": [{"type": "equidistant", "bins": 2}] * 2,
                        }
                    },
                }
            ],
        },
        "Volumes": {
            "acts-geometry-hierarchy-map": {
                "format-version": 0,
                "value-identifier": "Material Volume Map",
            },
            "entries": [],
        },
    }
    path = tmp_path / "prototype.json"
    path.write_text(json.dumps(document))
    acts.IMaterialDecorator.fromFile(path).decorate(surface)
    return surface


def test_grid_accumulator_bindings(tmp_path):
    surface = material_surface(tmp_path)
    cfg = acts.GridSurfaceMaterialAccumulator.Config(
        materialSurfaces=[surface], emptyBinCorrection=False
    )
    accumulator = acts.GridSurfaceMaterialAccumulator(cfg, acts.logging.INFO)
    assert isinstance(accumulator, acts.ISurfaceMaterialAccumulator)
    context = acts.GeometryContext.dangerouslyDefaultConstruct()
    state = accumulator.createState(context)
    assert isinstance(state, acts.ISurfaceMaterialAccumulator.State)
    accumulator.accumulate(state, context, [], [])
    assert len(accumulator.finalizeMaterial(state, context)) == 1
    maps = accumulator.finalizeMaps(state, context)
    output = tmp_path / "empty-grid.json"
    TrackingGeometryMaterialJsonConverter().toFile(maps, output)
    material = json.loads(output.read_text())["surfaces"][0]["material"]
    assert material["kind"] == "grid"
    assert material["storage"]["kind"] == "direct"
    assert len(material["storage"]["values"]) == 16
    assert [axis["bins"] for axis in material["axes"]] == [2, 2]


def write_material_tracks(path):
    ROOT = pytest.importorskip("ROOT")
    output = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree("material_tracks", "Deterministic material tracks")
    event = array("I", [0])
    tree.Branch("event_id", event, "event_id/i")
    values = {
        "v_x": -5.0,
        "v_y": -5.0,
        "v_z": -10.0,
        "v_px": 0.0,
        "v_py": 0.0,
        "v_pz": 1.0,
        "v_phi": 0.0,
        "v_eta": 0.0,
        "t_X0": 0.0,
        "t_L0": 0.0,
    }
    scalars = {name: array("f", [value]) for name, value in values.items()}
    for name, value in scalars.items():
        tree.Branch(name, value, name + "/F")
    steps = {
        "mat_x": -5.0,
        "mat_y": -5.0,
        "mat_z": 0.0,
        "mat_dx": 0.0,
        "mat_dy": 0.0,
        "mat_dz": 1.0,
        "mat_step_length": 0.0,
        "mat_X0": 10.0,
        "mat_L0": 20.0,
        "mat_A": 12.0,
        "mat_Z": 6.0,
        "mat_rho": 1.0,
    }
    vectors = {name: ROOT.std.vector("float")() for name in steps}
    for name, value in vectors.items():
        tree.Branch(name, value)
    elements = ROOT.std.vector("std::vector<unsigned int>")()
    fractions = ROOT.std.vector("std::vector<float>")()
    tree.Branch("elements", elements)
    tree.Branch("fraction", fractions)
    for index, thickness in enumerate((2.0, 6.0, 0.0)):
        event[0] = index
        scalars["t_X0"][0] = thickness / 10.0
        scalars["t_L0"][0] = thickness / 20.0
        for name, vector in vectors.items():
            vector.clear()
            if thickness:
                vector.push_back(
                    thickness if name == "mat_step_length" else steps[name]
                )
        tree.Fill()
    tree.Write()
    output.Close()


@pytest.mark.parametrize("accumulator", ["binned", "grid"])
def test_material_mapping_accumulator_switch(tmp_path, accumulator):
    from material_mapping import runMaterialMapping

    surface = material_surface(tmp_path)
    input_file = tmp_path / "input.root"
    write_material_tracks(input_file)
    output = tmp_path / accumulator
    sequencer = runMaterialMapping(
        [surface],
        input_file,
        str(output),
        loglevel=acts.logging.WARNING,
        accumulator=accumulator,
    )
    sequencer.run()
    # MaterialMapping writes the finalized map on destruction.
    del sequencer
    material = json.loads((tmp_path / f"{accumulator}_map.json").read_text())[
        "Surfaces"
    ]["entries"][0]["value"]["material"]
    assert material["type"] == ("grid" if accumulator == "grid" else "binned")
    if accumulator == "grid":
        data = material["accessor"]["grid"]["data"]
        slab = next(slab for bins, slab in data if bins == [1, 1])
        assert not (tmp_path / "grid_map.root").exists()
    else:
        slab = material["data"][0][0]
        assert (tmp_path / "binned_map.root").exists()
    assert slab["thickness"] == pytest.approx(8.0 / 3.0)
    assert (tmp_path / f"{accumulator}_mapped.root").exists()
    assert (tmp_path / f"{accumulator}_unmapped.root").exists()


def test_material_mapping_rejects_unsupported_choices(tmp_path):
    from material_mapping import runMaterialMapping

    with pytest.raises(ValueError, match="accumulator"):
        runMaterialMapping(
            [], tmp_path / "missing.root", "unused", accumulator="invalid"
        )
    with pytest.raises(ValueError, match="JSON or CBOR"):
        runMaterialMapping(
            [],
            tmp_path / "missing.root",
            "unused",
            outputMapFormats=["root"],
            accumulator="grid",
        )


@pytest.mark.parametrize(
    "accumulator,schema,suffix",
    [
        ("binned", "legacy", ".root"),
        ("binned", "legacy", ".json"),
        ("grid", "legacy", ".json"),
        ("grid", "legacy", ".cbor"),
        ("grid", "versioned", ".json"),
        ("grid", "versioned", ".cbor"),
        ("grid", "versioned", ".json.zst"),
        ("grid", "versioned", ".cbor.zst"),
    ],
)
def test_material_validation_reads_mapped_material(
    tmp_path, accumulator, schema, suffix
):
    from material_mapping import runMaterialMapping
    from material_validation import loadMaterialDecorator, runMaterialValidation

    ROOT = pytest.importorskip("ROOT")
    source = material_surface(tmp_path)
    input_file = tmp_path / "input.root"
    write_material_tracks(input_file)
    formats = ["json", "root"] if suffix == ".root" else ["json", "cbor"]
    sequencer = runMaterialMapping(
        [source],
        input_file,
        str(tmp_path / "mapped"),
        outputMapFormats=formats,
        loglevel=acts.logging.WARNING,
        accumulator=accumulator,
    )
    sequencer.run()
    del sequencer

    map_file = tmp_path / ("mapped_map" + suffix)
    if schema == "versioned":
        maps = loadMaterialDecorator(tmp_path / "mapped_map.json").materialMaps
        map_file = tmp_path / ("versioned_map" + suffix)
        try:
            TrackingGeometryMaterialJsonConverter().toFile(maps, map_file)
        except RuntimeError as error:
            if suffix.endswith(".zst") and "without zstd support" in str(error):
                pytest.skip("ACTS was built without zstd support")
            raise

    # The mapped bin is at local (-5, -5). Move the target plane forward so
    # fixed particle-gun directions from the origin cross that bin.
    target = acts.Surface.createPlane(
        acts.Transform3(acts.Vector3(0.0, 0.0, 10.0)), acts.RectangleBounds(10.0, 10.0)
    )
    target.assignGeometryId(source.geometryId)
    loadMaterialDecorator(map_file).decorate(target)
    assert target.surfaceMaterial is not None

    eta = math.asinh(math.sqrt(2.0))
    phi = 225.0 * acts.UnitConstants.degree
    output = tmp_path / "validated"
    sequencer = acts.examples.Sequencer(events=2, numThreads=1)
    runMaterialValidation(
        surfaces=[target],
        s=sequencer,
        tracksPerEvent=3,
        etaRange=(eta, eta),
        phiRange=(phi, phi),
        outputFileBase=output,
    ).run()
    del sequencer

    # The input tracks map to thickness 8/3; validation adds the incidence
    # correction sqrt(1 + 5^2/10^2 + 5^2/10^2).
    expected = (8.0 / 3.0) * math.sqrt(1.5)
    result = ROOT.TFile.Open(str(output) + ".root")
    tree = result.Get("material_tracks")
    assert tree.GetEntries() == 6
    for track in tree:
        assert list(track.mat_step_length) == pytest.approx([expected])
        assert track.t_X0 == pytest.approx(expected / 10.0)
        assert track.t_L0 == pytest.approx(expected / 20.0)
    result.Close()


def test_material_validation_rejects_unknown_map(tmp_path):
    from material_validation import loadMaterialDecorator

    path = tmp_path / "unknown.json"
    path.write_text(json.dumps({"surfaces": []}))
    with pytest.raises(ValueError, match="Unrecognized material map format"):
        loadMaterialDecorator(path)


def test_material_validation_rejects_unsupported_version(tmp_path):
    from material_validation import loadMaterialDecorator

    path = tmp_path / "unsupported.json"
    path.write_text(
        json.dumps(
            {
                "header": {"format": "acts-material-map", "version": 2},
                "surfaces": [],
            }
        )
    )
    with pytest.raises(ValueError, match="unsupported material document version"):
        loadMaterialDecorator(path)
