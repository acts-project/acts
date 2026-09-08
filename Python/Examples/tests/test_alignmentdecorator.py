import acts
import pytest

from acts.examples import Sequencer
from helpers import alignmentEnabled


@pytest.mark.skipif(not alignmentEnabled, reason="Alignment module is not set up")
def test_alignmentdecorator_io_mode(capfd):
    """This tests the alignment decorator in IO mode,
    i.e. with a given pre-defined alignment store"""

    from acts.examples.alignment import AlignmentDecorator, GeoIdAlignmentStore

    alignDecoConfig = AlignmentDecorator.Config()

    # Create a dummy alignment store
    geoId = acts.GeometryIdentifier(volume=1, layer=2)
    trf = acts.Transform3(acts.Vector3(0.0, 0.0, 0.0))
    geoIdMap = {}
    geoIdMap[geoId] = trf
    alignmentStore = GeoIdAlignmentStore(geoIdMap)
    alignDecoConfig.iovStores = [((10, 20), alignmentStore)]
    alignDecoConfig.garbageCollection = True
    alignDecoConfig.gcInterval = 20

    alignDeco = AlignmentDecorator(alignDecoConfig, acts.logging.DEBUG)

    sequencer = Sequencer(
        events=100,
        numThreads=4,
        logLevel=acts.logging.INFO,
    )

    sequencer.addContextDecorator(alignDeco)
    sequencer.run()
    if capfd is not None:
        out, err = capfd.readouterr()
        # Check that the alignment store is decorated for events 10 to 20
        for event in range(10, 20):
            assert (
                f"Decorating AlgorithmContext with alignment store for event {event}"
                in out
            )
        # Count that there is only one garbage collection call
        assert out.count("Garbage collection: removing alignment store") == 1


@pytest.mark.skipif(not alignmentEnabled, reason="Alignment module is not set up")
def test_alignmentdecorator_gen_mode(capfd):

    from acts.examples.alignment import (
        AlignmentDecorator,
        AlignmentGeneratorGlobalShift,
        AlignmentGeneratorGlobalRotation,
        GeoIdAlignmentStore,
    )

    """This tests the alignment decorator in generative mode"""
    alignDecoConfig = AlignmentDecorator.Config()

    # Create some nominal store
    geoId0 = acts.GeometryIdentifier(volume=1, layer=2)
    trf0 = acts.Transform3(acts.Vector3(0.0, 0.0, 10.0))
    geoId1 = acts.GeometryIdentifier(volume=1, layer=4)
    trf1 = acts.Transform3(acts.Vector3(0.0, 0.0, 20.0))
    geoIdMap = {}
    geoIdMap[geoId0] = trf0
    geoIdMap[geoId1] = trf1
    alignDecoConfig.nominalStore = GeoIdAlignmentStore(geoIdMap)

    gShift = AlignmentGeneratorGlobalShift()
    gShift.shift = acts.Vector3(0.0, 0.0, 100.0)

    gRot = AlignmentGeneratorGlobalRotation()
    gRot.axis = acts.Vector3(1.0, 0.0, 0.0)
    gRot.angle = 0.15

    alignDecoConfig.iovGenerators = [((10, 20), gShift), ((50, 75), gRot)]
    alignDecoConfig.garbageCollection = True
    alignDecoConfig.gcInterval = 20

    alignDeco = AlignmentDecorator(alignDecoConfig, acts.logging.VERBOSE)

    sequencer = Sequencer(
        events=100,
        numThreads=1,
        logLevel=acts.logging.INFO,
    )

    sequencer.addContextDecorator(alignDeco)
    sequencer.run()
    # Count that the alignment store is decorated 37 times
    out, err = capfd.readouterr()
    assert out.count("Decorating AlgorithmContext with alignment store") == 37


def _meanResidualPerLayer(rootFile, layer):
    """Mean smoothed local residual on one telescope layer, and on all others."""
    import uproot
    import numpy as np

    tree = uproot.open(f"{rootFile}:trackstates")
    data = tree.arrays(["volume_id", "layer_id", "res_eLOC1_smt"], library="np")
    vol = np.concatenate(data["volume_id"])
    lay = np.concatenate(data["layer_id"])
    res = np.concatenate(data["res_eLOC1_smt"])

    onLayer = (vol == 1) & (lay == layer) & np.isfinite(res)
    elsewhere = (vol == 1) & (lay != layer) & np.isfinite(res)
    assert onLayer.sum() > 50, "not enough track states on the misaligned layer"
    return np.mean(res[onLayer]), np.mean(np.abs(res[elsewhere]))


@pytest.mark.skipif(not alignmentEnabled, reason="Alignment module is not set up")
@pytest.mark.parametrize("target", ["both", "sim"])
def test_alignmentdecorator_sim_reco_split(tmp_path, target):
    """Simulation and reconstruction can run on different module placements.

    With target=eBoth the shift is known to reconstruction and residuals stay
    centred. With target=eSim only simulation sees it and the shifted layer picks
    up a bias.
    """
    pytest.importorskip("uproot")

    from acts import UnitConstants as u
    from acts.examples import RandomNumbers, Sequencer, TelescopeDetector
    from acts.examples.alignment import AlignmentDecorator

    from misaligned_simulation import (
        LAYER_POSITIONS,
        addMisalignment,
        addTelescopeChain,
    )

    layer = 4
    shift = 0.2 * u.mm

    detector = TelescopeDetector(
        bounds=[200, 200],
        positions=LAYER_POSITIONS,
        stereos=[0] * len(LAYER_POSITIONS),
        binValue=1,
    )
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    s = Sequencer(events=300, numThreads=1, outputDir=str(tmp_path))
    addMisalignment(
        s,
        trackingGeometry,
        layer=layer,
        shift=acts.Vector3(0, 0, shift),
        target={
            "sim": AlignmentDecorator.Target.eSim,
            "both": AlignmentDecorator.Target.eBoth,
        }[target],
    )
    addTelescopeChain(s, trackingGeometry, field, tmp_path, RandomNumbers(seed=42))
    s.run()

    onLayer, elsewhere = _meanResidualPerLayer(tmp_path / "trackstates_ckf.root", layer)

    if target == "both":
        # reconstruction knows where the layer is, nothing is biased
        assert abs(onLayer) < 0.25 * shift
    else:
        # the fit absorbs part of the shift by tilting the track, so the bias on
        # the shifted layer is a sizeable fraction of it rather than all of it
        assert abs(onLayer) > 0.5 * shift
        assert abs(onLayer) > 2 * elsewhere
