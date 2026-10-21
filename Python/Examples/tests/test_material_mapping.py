from pathlib import Path
import json
import os
import shutil
import subprocess


from helpers import (
    geant4Enabled,
    dd4hepEnabled,
    assert_entries,
)
import pytest

import acts
from acts.examples.odd import getOpenDataDetector

from acts.examples import Sequencer


@pytest.mark.slow
@pytest.mark.odd
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_material_recording(tmp_path, material_recording, assert_root_hash):
    root_files = [
        (
            "geant4_material_tracks.root",
            "material_tracks",
            2000,
        )
    ]

    for fn, tn, ee in root_files:
        fp = material_recording / fn
        assert fp.exists()
        assert fp.stat().st_size > 2**10 * 50
        assert_entries(fp, tn, ee)
        assert_root_hash(fn, fp)


@pytest.mark.slow
@pytest.mark.odd
@pytest.mark.parametrize(
    "threads,region",
    [(1, None), (2, None), (1, (500.0, -300.0, 800.0))],
    ids=["st", "mt", "region"],
)
def test_material_recording_ddsim(tmp_path, threads, region):
    if shutil.which("ddsim") is None:
        pytest.skip("ddsim not available")
    # ddsim finds the DDG4 actions through the .components file of their
    # library on the library path
    if not any(
        Path(d, "libActsPluginDD4hepGeant4.components").exists()
        for d in os.environ.get("LD_LIBRARY_PATH", "").split(os.pathsep)
        if d
    ):
        pytest.skip("DDG4 material recording actions not built")
    if threads > 1:
        # ddsim only runs multi-threaded from DD4hep 1.38
        ddsimHelp = subprocess.run(
            ["ddsim", "--help"], capture_output=True, text=True, check=False
        )
        if "--numberOfThreads" not in ddsimHelp.stdout:
            pytest.skip("ddsim does not support multi-threading")

    from acts.examples.odd import getOpenDataDetectorDirectory
    from acts.examples.root import RootMaterialTrackReader, RootMaterialTrackWriter

    steering = (
        Path(__file__).parent.parent.parent.parent
        / "Examples/Scripts/Python/material_recording_ddsim.py"
    )
    oddXml = getOpenDataDetectorDirectory() / "xml"
    output = tmp_path / "ddsim_material_tracks.root"
    events = 4
    tracks = 50

    # The tracker alone keeps the geometry construction fast
    cmd = [
        "ddsim",
        "--steeringFile",
        str(steering),
        "--compactFile",
        str(oddXml / "OpenDataDetectorDefs.xml"),
        str(oddXml / "OpenDataDetectorTracker.xml"),
        "--numberOfEvents",
        str(events),
        "--gun.multiplicity",
        str(tracks),
        "--random.seed",
        "42",
        "--outputFile",
        str(output),
    ]
    if threads > 1:
        cmd += ["--numberOfThreads", str(threads)]
    if region is not None:
        # An asymmetric region in z
        maxR, minZ, maxZ = region
        cmd += [
            f"--outputConfig.recordingMaxR={maxR}",
            f"--outputConfig.recordingMinZ={minZ}",
            f"--outputConfig.recordingMaxZ={maxZ}",
        ]
    result = subprocess.run(
        cmd, cwd=tmp_path, capture_output=True, text=True, check=False
    )
    assert (
        result.returncode == 0
    ), f"ddsim failed:\n{result.stdout[-5000:]}\n{result.stderr[-5000:]}"

    # One entry per geantino, also for those that cross no material
    assert output.exists()
    assert_entries(output, "material_tracks", events * tracks)
    import ROOT

    rdf = ROOT.RDataFrame("material_tracks", str(output))
    data = rdf.AsNumpy(["event_id", "t_X0"])
    assert set(data["event_id"]) == set(range(events))
    if region is None:
        # Without a region, every geantino crosses at least the beam pipe
        assert (data["t_X0"] > 0).all()
    else:
        # Forward geantinos can leave the region before the beam pipe
        assert (data["t_X0"] > 0).mean() > 0.5

    if region is not None:
        # No recorded step ends outside the recording region
        ends = rdf.Define(
            "endR", "ROOT::VecOps::sqrt(mat_ex * mat_ex + mat_ey * mat_ey)"
        ).AsNumpy(["endR", "mat_ez"])
        steps = [i for i, r in enumerate(ends["endR"]) if len(r) > 0]
        tolerance = 1e-3
        assert max(ends["endR"][i].max() for i in steps) <= maxR + tolerance
        assert min(ends["mat_ez"][i].min() for i in steps) >= minZ - tolerance
        assert max(ends["mat_ez"][i].max() for i in steps) <= maxZ + tolerance

    # Read back with the reader the material mapping uses
    s = Sequencer(numThreads=1)
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks="material_tracks",
            treeName="material_tracks",
            fileList=[str(output)],
            readCachedSurfaceInformation=False,
        )
    )
    roundTrip = tmp_path / "round_trip.root"
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks="material_tracks",
            treeName="material_tracks",
            filePath=str(roundTrip),
            prePostStep=True,
            recalculateTotals=True,
        )
    )
    s.run()
    assert_entries(roundTrip, "material_tracks", events * tracks)


@pytest.mark.slow
@pytest.mark.odd
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.parametrize("gen3", [False, True], ids=["gen1", "gen3"])
def test_material_mapping(
    material_recording, tmp_path, assert_root_hash, request, gen3
):
    from material_mapping import runMaterialMapping
    from material_validation import runMaterialValidation

    map_file = tmp_path / "material-map_tracks.root"
    assert not map_file.exists()

    odd_detector = request.getfixturevalue(
        "odd_detector_gen3" if gen3 else "odd_detector_gen1"
    )
    trackingGeometry = odd_detector.trackingGeometry()
    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    s = Sequencer(events=2000, numThreads=1)

    s, mapping = runMaterialMapping(
        surfaces=materialSurfaces,
        inputFile=material_recording / "geant4_material_tracks.root",
        outputFileBase=str(tmp_path / "material_mapping"),
        loglevel=acts.logging.INFO,
        outputMaterialTracks="material_tracks",
        treeName="material_tracks",
        s=s,
    )

    s.run()

    from acts.json import TrackingGeometryMaterialJsonConverter
    from acts.examples.root import RootMaterialWriter

    TrackingGeometryMaterialJsonConverter().toFile(
        mapping.material, tmp_path / "material_mapping_map.json"
    )
    # Gen3 faces carry stable material keys, which the ROOT material-map
    # writer rejects, so Gen3 maps are JSON-only.
    if not gen3:
        RootMaterialWriter(
            level=acts.logging.INFO,
            filePath=str(tmp_path / "material_mapping_map.root"),
        ).writeMaterial(mapping.material)

        # root map output check
        map_file_root = tmp_path / "material_mapping_map.root"
        assert map_file_root.exists()
        assert_root_hash(map_file_root.name, map_file_root)

    # json map output check
    map_file_json = tmp_path / "material_mapping_map.json"
    assert map_file_json.exists()
    with map_file_json.open() as fh:
        assert json.load(fh)["header"]["version"] == 1

    # mapped tracks output check
    map_file_mapped = tmp_path / "material_mapping_mapped.root"
    assert map_file_mapped.exists()
    assert_entries(map_file_mapped, "material_tracks", 2000)
    assert_root_hash(map_file_mapped.name, map_file_mapped)

    # unmapped tracks output check
    map_file_unmapped = tmp_path / "material_mapping_unmapped.root"
    assert map_file_unmapped.exists()
    assert_entries(map_file_unmapped, "material_tracks", 2000)
    assert_root_hash(map_file_unmapped.name, map_file_unmapped)

    val_file = tmp_path / "material_validation.root"
    assert not val_file.exists()

    # test the validation as well
    s = Sequencer(events=10, numThreads=1)

    with getOpenDataDetector(
        materialDecorator=acts.IMaterialDecorator.fromFile(map_file_json), gen3=gen3
    ) as detector:
        trackingGeometry = detector.trackingGeometry()
        materialSurfaces = trackingGeometry.extractMaterialSurfaces()

        runMaterialValidation(
            surfaces=materialSurfaces,
            s=s,
            tracksPerEvent=1000,
            outputFileBase=tmp_path / "material_validation",
            materialTrackCollectionName="material_tracks",
        )

        s.run()

    assert val_file.exists()
    assert_entries(val_file, "material_tracks", 10000)
    assert_root_hash(val_file.name, val_file)
