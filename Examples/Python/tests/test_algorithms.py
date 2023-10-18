import pytest

import acts

from acts.examples import (
    TutorialVertexFinderAlgorithm,
    AdaptiveMultiVertexFinderAlgorithm,
    VertexFitterAlgorithm,
    IterativeVertexFinderAlgorithm,
    SpacePointMaker,
    TrackFindingAlgorithm,
    SeedingAlgorithm,
    TrackParamsEstimationAlgorithm,
    EventGenerator,
    FatrasSimulation,
    MaterialMapping,
    TruthSeedSelector,
    TruthTrackFinder,
    ParticleSelector,
    TruthVertexFinder,
    ParticleSmearing,
    TrackSelectorAlgorithm,
    TrackFittingAlgorithm,
    SurfaceSortingAlgorithm,
    ParticlesPrinter,
    HitsPrinter,
    TrackParametersPrinter,
    PropagationAlgorithm,
    DigitizationAlgorithm,
    PlanarSteppingAlgorithm,
)


from helpers import geant4Enabled, hepmc3Enabled


@pytest.mark.parametrize(
    "alg",
    [
        TutorialVertexFinderAlgorithm,
        AdaptiveMultiVertexFinderAlgorithm,
        VertexFitterAlgorithm,
        IterativeVertexFinderAlgorithm,
        SpacePointMaker,
        TrackFindingAlgorithm,
        SeedingAlgorithm,
        TrackParamsEstimationAlgorithm,
        EventGenerator,
        FatrasSimulation,
        MaterialMapping,
        TruthSeedSelector,
        TruthTrackFinder,
        ParticleSelector,
        TruthVertexFinder,
        ParticleSmearing,
        TrackSelectorAlgorithm,
        TrackFittingAlgorithm,
        SurfaceSortingAlgorithm,
        ParticlesPrinter,
        HitsPrinter,
        TrackParametersPrinter,
        PropagationAlgorithm,
        # GeantinoRecording,
        PlanarSteppingAlgorithm,
        # EventRecording,
    ],
)
def test_algorithm_interface(alg):
    assert hasattr(alg, "Config")


@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
@pytest.mark.skipif(not hepmc3Enabled, reason="HepMC3 not set up")
def test_g4_algorithms():
    from acts.examples.geant4.hepmc3 import EventRecording
    from acts.examples.geant4 import Geant4Simulation

    for alg in (EventRecording, Geant4Simulation):
        assert hasattr(alg, "Config")


@pytest.mark.skipif(not hepmc3Enabled, reason="HepMC3 not set up")
def test_hepmc_algorithms():
    from acts.examples.hepmc3 import HepMCProcessExtractor

    assert hasattr(HepMCProcessExtractor, "Config")


def test_special_algorithm_interfaces():
    # just assert they exists
    assert DigitizationAlgorithm
