import pytest

from acts.examples import (
    AdaptiveMultiVertexFinderAlgorithm,
    VertexFitterAlgorithm,
    IterativeVertexFinderAlgorithm,
    HoughVertexFinderAlgorithm,
    SpacePointMaker,
    TrackFindingAlgorithm,
    SeedingAlgorithm,
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
        SeedingAlgorithm,
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
