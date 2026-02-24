import pytest

import acts
from acts.examples import (
    AdaptiveMultiVertexFinderAlgorithm,
    VertexFitterAlgorithm,
    IterativeVertexFinderAlgorithm,
    HoughVertexFinderAlgorithm,
    SpacePointMaker,
    TrackFindingAlgorithm,
    GridTripletSeedingAlgorithm,
    OrthogonalTripletSeedingAlgorithm,
    TrackParamsEstimationAlgorithm,
    EventGenerator,
    FatrasSimulation,
    MaterialMapping,
    TruthTrackFinder,
    ParticleSelector,
    TruthVertexFinder,
    TrackParameterSmearing,
    TrackSelectorAlgorithm,
    TrackFittingAlgorithm,
    ParticlesPrinter,
    TrackParametersPrinter,
    PropagationAlgorithm,
    DigitizationAlgorithm,
)


from helpers import geant4Enabled, hepmc3Enabled


@pytest.mark.parametrize(
    "alg",
    [
        AdaptiveMultiVertexFinderAlgorithm,
        VertexFitterAlgorithm,
        IterativeVertexFinderAlgorithm,
        HoughVertexFinderAlgorithm,
        SpacePointMaker,
        TrackFindingAlgorithm,
        GridTripletSeedingAlgorithm,
        OrthogonalTripletSeedingAlgorithm,
        TrackParamsEstimationAlgorithm,
        EventGenerator,
        FatrasSimulation,
        MaterialMapping,
        TruthTrackFinder,
        ParticleSelector,
        TruthVertexFinder,
        TrackParameterSmearing,
        TrackSelectorAlgorithm,
        TrackFittingAlgorithm,
        ParticlesPrinter,
        TrackParametersPrinter,
        PropagationAlgorithm,
        # GeantinoRecording,
        # EventRecording,
    ],
)
def test_algorithm_interface(alg):
    assert hasattr(alg, "Config")


@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_g4_algorithms():
    from acts.examples.geant4 import Geant4Simulation

    assert hasattr(Geant4Simulation, "Config")


def test_special_algorithm_interfaces():
    # just assert they exists
    assert DigitizationAlgorithm


def test_ialgorithm_init_with_log_level():
    """PyIAlgorithm subclass can be initialized with name and log level only."""

    class AlgWithLevel(acts.examples.IAlgorithm):
        def __init__(self):
            super().__init__(name="AlgWithLevel", level=acts.logging.DEBUG)

        def execute(self, context):
            return acts.examples.ProcessCode.SUCCESS

    alg = AlgWithLevel()
    assert alg.name() == "AlgWithLevel"
    assert alg.logger.level == acts.logging.DEBUG


def test_ialgorithm_init_with_logger():
    """PyIAlgorithm subclass can be initialized with name and a full logger object."""

    class AlgWithLogger(acts.examples.IAlgorithm):
        def __init__(self):
            logger = acts.getDefaultLogger("AlgWithLogger", acts.logging.VERBOSE)
            super().__init__(name="AlgWithLogger", logger=logger)

        def execute(self, context):
            return acts.examples.ProcessCode.SUCCESS

    alg = AlgWithLogger()
    assert alg.name() == "AlgWithLogger"
    assert alg.logger.level == acts.logging.VERBOSE
    assert alg.logger.name == "AlgWithLogger"
