"""Map and validate binned/grid material with deterministic material tracks."""

from array import array
import json
import math
import shutil
import sys
from types import SimpleNamespace

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
@pytest.mark.parametrize("custom_sequencer", [False, True])
def test_material_mapping_accumulator_switch(tmp_path, accumulator, custom_sequencer):
    from material_mapping import runMaterialMapping

    surface = material_surface(tmp_path)
    input_file = tmp_path / "input.root"
    write_material_tracks(input_file)
    output = tmp_path / accumulator
    supplied_sequencer = (
        acts.examples.Sequencer(events=2, numThreads=1) if custom_sequencer else None
    )
    sequencer, mapping = runMaterialMapping(
        [surface],
        input_file,
        str(output),
        loglevel=acts.logging.WARNING,
        s=supplied_sequencer,
        accumulator=accumulator,
    )
    if custom_sequencer:
        assert sequencer is supplied_sequencer
    sequencer.run()
    TrackingGeometryMaterialJsonConverter().toFile(
        mapping.material, tmp_path / f"{accumulator}_map.json"
    )
    material = json.loads((tmp_path / f"{accumulator}_map.json").read_text())[
        "surfaces"
    ][0]["material"]
    assert material["kind"] == accumulator
    if accumulator == "grid":
        slabs = material["storage"]["values"]
        assert not (tmp_path / "grid_map.root").exists()
    else:
        slabs = material["values"]
        from acts.examples.root import RootMaterialWriter

        RootMaterialWriter(
            level=acts.logging.WARNING, filePath=str(output) + "_map.root"
        ).writeMaterial(mapping.material)
        assert (tmp_path / "binned_map.root").exists()
    assert max(slab["thickness"] for slab in slabs) == pytest.approx(
        4.0 if custom_sequencer else 8.0 / 3.0
    )
    assert (tmp_path / f"{accumulator}_mapped.root").exists()
    assert (tmp_path / f"{accumulator}_unmapped.root").exists()


def test_material_mapping_rejects_unsupported_choices(tmp_path, monkeypatch, capsys):
    import material_mapping
    from material_mapping import runMaterialMapping

    with pytest.raises(ValueError, match="accumulator"):
        runMaterialMapping(
            [], tmp_path / "missing.root", "unused", accumulator="invalid"
        )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "material_mapping.py",
            "--accumulator",
            "grid",
            "--output-map-formats",
            "root",
        ],
    )
    with pytest.raises(SystemExit) as error:
        material_mapping.main()
    assert error.value.code == 2
    assert "JSON or CBOR" in capsys.readouterr().err


@pytest.mark.parametrize("explicit_config", [False, True])
@pytest.mark.parametrize("accumulator", ["binned", "grid"])
@pytest.mark.parametrize("formats", [None, ["json", "cbor"]])
def test_mapping_cli_material_configuration(
    tmp_path, monkeypatch, explicit_config, accumulator, formats
):
    import material_mapping
    from material_validation import loadMaterialDecorator

    surface = material_surface(tmp_path)
    odd_config = tmp_path / "config"
    odd_config.mkdir()
    shutil.copy(
        tmp_path / "prototype.json", odd_config / "odd-material-mapping-config.json"
    )
    input_file = tmp_path / "input.root"
    write_material_tracks(input_file)

    # Exercise actual prototype loading and mapping at the CLI boundary while
    # replacing the expensive ODD construction with one material surface.
    def detector(materialDecorator):
        if explicit_config or accumulator == "grid":
            assert materialDecorator is not None
            materialDecorator.decorate(surface)
        else:
            assert materialDecorator is None
        geometry = SimpleNamespace(extractMaterialSurfaces=lambda: [surface])
        return SimpleNamespace(trackingGeometry=lambda: geometry)

    def odd_directory():
        assert not explicit_config, "Explicit --matconfig must take precedence"
        assert accumulator == "grid", "Only grid mode defaults to prototype material"
        return tmp_path

    monkeypatch.setattr(material_mapping, "getOpenDataDetector", detector)
    monkeypatch.setattr(material_mapping, "getOpenDataDetectorDirectory", odd_directory)
    arguments = [
        "material_mapping.py",
        "--accumulator",
        accumulator,
        "-i",
        str(input_file),
        "-o",
        str(tmp_path / "mapped"),
        "-n2",
    ]
    if formats is not None:
        arguments += ["--output-map-formats", *formats]
    if explicit_config:
        arguments += ["--matconfig", str(tmp_path / "prototype.json")]
    monkeypatch.setattr(sys, "argv", arguments)
    material_mapping.main()

    suffix = "json" if formats is None else "cbor"
    maps = loadMaterialDecorator(tmp_path / f"mapped_map.{suffix}").materialMaps
    json_map = tmp_path / "mapped-versioned.json"
    TrackingGeometryMaterialJsonConverter().toFile(maps, json_map)
    material = json.loads(json_map.read_text())["surfaces"][0]["material"]
    assert material["kind"] == accumulator
    # -n2 processes only the two non-empty input events; including the third
    # empty crossing would change the mapped thickness from 4 to 8/3.
    slabs = (
        material["storage"]["values"] if accumulator == "grid" else material["values"]
    )
    assert max(slab["thickness"] for slab in slabs) == 4.0
    assert (tmp_path / "mapped_map.root").exists() == (
        accumulator == "binned" and formats is None
    )
    assert (tmp_path / "mapped_map.cbor").exists() == (formats is not None)


@pytest.mark.parametrize(
    "accumulator,schema,suffix",
    [
        ("binned", "legacy", ".root"),
        ("binned", "legacy", ".json"),
        ("binned", "legacy", ".cbor"),
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
    sequencer, mapping = runMaterialMapping(
        [source],
        input_file,
        str(tmp_path / "mapped"),
        loglevel=acts.logging.WARNING,
        accumulator=accumulator,
    )
    sequencer.run()

    map_file = tmp_path / ("mapped_map" + suffix)
    if suffix == ".root":
        from acts.examples.root import RootMaterialWriter

        RootMaterialWriter(
            level=acts.logging.WARNING, filePath=str(map_file)
        ).writeMaterial(mapping.material)
    elif schema == "legacy":
        from acts.examples.json import JsonFormat, JsonMaterialWriter
        from acts.json import MaterialMapJsonConverter

        JsonMaterialWriter(
            level=acts.logging.WARNING,
            converterCfg=MaterialMapJsonConverter.Config(
                processSensitives=True,
                processApproaches=True,
                processRepresenting=True,
                processBoundaries=True,
                processVolumes=False,
            ),
            fileName=str(tmp_path / "mapped_map"),
            writeFormat=JsonFormat.Cbor if suffix == ".cbor" else JsonFormat.Json,
        ).writeMaterial(mapping.material)
    else:
        try:
            TrackingGeometryMaterialJsonConverter().toFile(mapping.material, map_file)
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
