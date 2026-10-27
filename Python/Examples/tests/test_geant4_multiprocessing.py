from functools import partial
import json
import multiprocessing
import os

import pytest

from helpers import dd4hepEnabled, geant4Enabled, hepmc3Enabled, rootEnabled
from acts.examples.geant4_multiprocessing import runGeant4Processes


def _write_event_range(begin, end, *, directory, barrier):
    # All workers must reach the barrier together: sequential subprocesses
    # cannot pass this check.
    barrier.wait(timeout=30)
    (directory / f"{begin}.json").write_text(
        json.dumps(
            dict(
                pid=os.getpid(),
                startMethod=multiprocessing.get_start_method(),
                events=list(range(begin, end)),
            )
        )
    )


def _fail_event_range(begin, end, *, hardExit):
    if begin == 0:
        if hardExit:
            os._exit(17)
        raise RuntimeError("worker failure")
    # Failure must terminate this worker rather than waiting for it to finish.
    multiprocessing.get_context("spawn").Event().wait(30)


def _telescope_detector():
    import acts.examples

    return acts.examples.TelescopeDetector(
        acts.examples.TelescopeDetector.Config(
            bounds=[200, 200], positions=[30, 60, 90], stereos=[0, 0, 0]
        )
    )


def _broken_detector():
    raise RuntimeError("detector construction failed")


def test_geant4_processes_run_concurrently(tmp_path):
    barrier = multiprocessing.get_context("spawn").Barrier(3)
    ranges = runGeant4Processes(
        partial(_write_event_range, directory=tmp_path, barrier=barrier),
        events=5,
        processes=3,
        skip=7,
    )
    records = [
        json.loads((tmp_path / f"{begin}.json").read_text()) for begin, _ in ranges
    ]
    assert len({r["pid"] for r in records}) == 3
    assert all(r["pid"] != os.getpid() and r["startMethod"] == "spawn" for r in records)
    assert [event for r in records for event in r["events"]] == list(range(7, 12))


@pytest.mark.parametrize("hardExit, exitCode", [(False, 1), (True, 17)])
def test_geant4_worker_failure_cleans_up_children(hardExit, exitCode):
    children = {child.pid for child in multiprocessing.active_children()}
    with pytest.raises(RuntimeError, match=f"failed with exit code {exitCode}"):
        runGeant4Processes(
            partial(_fail_event_range, hardExit=hardExit), events=2, processes=2
        )
    assert {child.pid for child in multiprocessing.active_children()} <= children


@pytest.mark.parametrize("options", [dict(events=-1), dict(processes=0), dict(skip=-1)])
def test_geant4_processes_reject_invalid_ranges(options):
    with pytest.raises(ValueError):
        runGeant4Processes(_fail_event_range, **(dict(events=1, processes=1) | options))


def _root_arrays(path):
    import uproot

    with uproot.open(path) as file:
        return {key: file[key].arrays(library="ak") for key in file.keys()}


@pytest.mark.skipif(
    not (geant4Enabled and hepmc3Enabled and rootEnabled and dd4hepEnabled),
    reason="Geant4 material recording dependencies not available",
)
@pytest.mark.parametrize("events, processes", [(5, 2), (1, 3), (0, 3)])
def test_parallel_material_recording_matches_serial(tmp_path, events, processes):
    import awkward as ak
    from material_recording import runMaterialRecordingParallel

    outputs = []
    for name, count in (("serial", 1), ("parallel", processes)):
        output = tmp_path / name
        runMaterialRecordingParallel(
            _telescope_detector,
            events=events,
            processes=count,
            skip=7,
            tracksPerEvent=10,
            outputFileBase=output,
        )
        outputs.append(_root_arrays(output.with_suffix(".root")))
    assert outputs[0].keys() == outputs[1].keys()
    for key in outputs[0]:
        assert len(outputs[0][key]) == events * 10
        assert ak.almost_equal(outputs[0][key], outputs[1][key], rtol=0, atol=0)


@pytest.mark.skipif(
    not (geant4Enabled and hepmc3Enabled and rootEnabled and dd4hepEnabled),
    reason="Geant4 material recording dependencies not available",
)
def test_failed_material_worker_preserves_output(tmp_path):
    from material_recording import runMaterialRecordingParallel

    output = tmp_path / "material.root"
    output.write_bytes(b"previous output")
    with pytest.raises(RuntimeError, match="failed with exit code"):
        runMaterialRecordingParallel(
            _broken_detector,
            events=3,
            processes=2,
            outputFileBase=output.with_suffix(""),
        )
    assert output.read_bytes() == b"previous output"


@pytest.mark.skipif(
    not (geant4Enabled and hepmc3Enabled and rootEnabled),
    reason="Geant4 simulation dependencies not available",
)
def test_parallel_geant4_matches_serial(tmp_path):
    import awkward as ak
    from geant4_parallel import runGeant4Parallel

    for name, count in (("serial", 1), ("parallel", 2)):
        runGeant4Parallel(
            events=5,
            processes=count,
            skip=7,
            outputDir=tmp_path / name,
            outputRoot=True,
            detectorFactory=_telescope_detector,
            etaRange=(4.0, 5.0),
        )
    serial = tmp_path / "serial"
    parallel = tmp_path / "parallel"
    files = sorted((serial / "csv").glob("*.csv"))
    assert files
    assert {p.name for p in files} == {p.name for p in (parallel / "csv").glob("*.csv")}
    for file in files:
        assert file.read_bytes() == (parallel / "csv" / file.name).read_bytes()
    for name in ("particles_simulation.root", "hits.root"):
        a, b = _root_arrays(serial / name), _root_arrays(parallel / name)
        assert a.keys() == b.keys()
        for key in a:
            assert len(a[key]) > 0
            assert ak.almost_equal(a[key], b[key], rtol=0, atol=0)
