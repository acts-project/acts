"""Check both accumulator choices with deterministic recorded material tracks."""

from array import array
import json

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
