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
