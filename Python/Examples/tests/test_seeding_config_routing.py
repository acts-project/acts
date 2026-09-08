# Regression tests for the field-forwarding bugs described in the "ACTS
# seeding configuration issues" writeup: addGridTripletSeeding /
# addOrthogonalTripletSeeding used to silently drop several config-arg
# fields instead of forwarding them to the algorithm's Config. Each test
# here sets every field the function is expected to forward to a
# distinctive, non-default value and asserts it lands on
# seedingAlg.config unchanged -- a silent drop now fails loudly instead of
# only showing up as two configs producing identical output.

import pytest

import acts
import acts.examples
from acts import UnitConstants as u
from acts.examples.reconstruction import (
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    addGridTripletSeeding,
    addOrthogonalTripletSeeding,
)


class _CapturingSequencer:
    """Stands in for acts.examples.Sequencer.

    addGridTripletSeeding/addOrthogonalTripletSeeding only ever call
    sequence.addAlgorithm(seedingAlg) and, via defaultLogging(),
    sequence.config.logLevel; capturing addAlgorithm's argument is the
    only way to get at the constructed algorithm's config, since
    Sequencer does not expose the algorithms it holds back to Python.
    """

    class _Config:
        logLevel = acts.logging.INFO

    def __init__(self):
        self.algorithm = None
        self.config = self._Config()

    def addAlgorithm(self, algorithm):
        self.algorithm = algorithm


def _confirmationRange(offset):
    return acts.SeedConfirmationRangeConfig(
        zMinSeedConf=-500 * u.mm - offset,
        zMaxSeedConf=500 * u.mm + offset,
        rMaxSeedConf=140 * u.mm + offset,
        nTopForLargeR=1,
        nTopForSmallR=2,
        seedConfMinBottomRadius=60 * u.mm + offset,
        seedConfMaxZOrigin=150 * u.mm + offset,
        minImpactSeedConf=1 * u.mm + offset,
    )


def _seedFinderConfigArg(**overrides):
    args = dict(
        maxSeedsPerSpM=7,
        cotThetaMax=1.11,
        sigmaScattering=2.22,
        radLengthPerSeed=0.033,
        minPt=0.44 * u.GeV,
        impactMax=1.55 * u.mm,
        interactionPointCut=True,
        deltaZMin=-77 * u.mm,
        deltaZMax=88 * u.mm,
        # zBinsCustomLooping is left at its default (empty): GridTriplet's
        # constructor requires every index in it to be < len(zBinEdges),
        # and zBinEdges is unset in most of these fixtures; the
        # forwarding test below sets both together explicitly.
        rRangeMiddleSP=[[10.0, 20.0]],
        useVariableMiddleSPRange=True,
        seedConfirmation=True,
        centralSeedConfirmationRange=_confirmationRange(0),
        forwardSeedConfirmationRange=_confirmationRange(1),
        rMinMiddle=12.3 * u.mm,
        rMaxMiddle=45.6 * u.mm,
        helixCutTolerance=2.5,
        toleranceParam=0.9 * u.mm,
        deltaR=(1 * u.mm, 2 * u.mm),
        deltaRBottomSP=(3 * u.mm, 4 * u.mm),
        deltaRTopSP=(5 * u.mm, 6 * u.mm),
        deltaRMiddleSPRange=(7 * u.mm, 8 * u.mm),
        collisionRegion=(-9 * u.mm, 9 * u.mm),
        r=(10 * u.mm, 20 * u.mm),
        z=(-30 * u.mm, 30 * u.mm),
    )
    args.update(overrides)
    return SeedFinderConfigArg(**args)


def _seedFilterConfigArg(**overrides):
    args = dict(
        impactWeightFactor=1.1,
        zOriginWeightFactor=2.2,
        compatSeedWeight=3.3,
        compatSeedLimit=4,
        numSeedIncrement=5.5,
        seedWeightIncrement=6.6,
        maxSeedsPerSpMConf=7,
        maxQualitySeedsPerSpMConf=8,
        useDeltaRorTopRadius=True,
        deltaInvHelixDiameter=0.0001 / u.mm,
    )
    args.update(overrides)
    return SeedFilterConfigArg(**args)


def _seedingAlgorithmConfigArg(**overrides):
    args = dict(
        zBinNeighborsTop=[(1, 2), (3, 4)],
        zBinNeighborsBottom=[(5, 6), (7, 8)],
        numPhiNeighbors=3,
        useExtraCuts=True,
    )
    args.update(overrides)
    return SeedingAlgorithmConfigArg(**args)


# Config fields are single-precision (C++ float); pytest.approx absorbs the
# double -> float rounding that a plain == would occasionally trip on.
def _approx(x):
    return pytest.approx(x, rel=1e-6)


# Fields common to both GridTripletSeedingAlgorithm::Config and
# OrthogonalTripletSeedingAlgorithm::Config, and forwarded by both
# functions the same way.
def _assertCommonFields(cfg, sfa, sfia, sfila):
    assert cfg.bFieldInZ == _approx(sfia.bFieldInZ)
    assert cfg.minPt == _approx(sfa.minPt)
    assert cfg.cotThetaMax == _approx(sfa.cotThetaMax)
    assert cfg.impactMax == _approx(sfa.impactMax)
    assert cfg.deltaRMin == _approx(sfa.deltaR[0])
    assert cfg.deltaRMax == _approx(sfa.deltaR[1])
    assert cfg.deltaRMinTop == _approx(sfa.deltaRTopSP[0])
    assert cfg.deltaRMaxTop == _approx(sfa.deltaRTopSP[1])
    assert cfg.deltaRMinBottom == _approx(sfa.deltaRBottomSP[0])
    assert cfg.deltaRMaxBottom == _approx(sfa.deltaRBottomSP[1])
    assert cfg.rMin == _approx(sfa.r[0])
    assert cfg.rMax == _approx(sfa.r[1])
    assert cfg.zMin == _approx(sfa.z[0])
    assert cfg.zMax == _approx(sfa.z[1])
    assert cfg.rMinMiddle == _approx(sfa.rMinMiddle)
    assert cfg.rMaxMiddle == _approx(sfa.rMaxMiddle)
    assert cfg.useVariableMiddleSPRange == sfa.useVariableMiddleSPRange
    assert list(cfg.rRangeMiddleSP) == sfa.rRangeMiddleSP
    assert cfg.deltaRMiddleMinSPRange == _approx(sfa.deltaRMiddleSPRange[0])
    assert cfg.deltaRMiddleMaxSPRange == _approx(sfa.deltaRMiddleSPRange[1])
    assert cfg.deltaZMin == _approx(sfa.deltaZMin)
    assert cfg.deltaZMax == _approx(sfa.deltaZMax)
    assert cfg.interactionPointCut == sfa.interactionPointCut
    assert cfg.collisionRegionMin == _approx(sfa.collisionRegion[0])
    assert cfg.collisionRegionMax == _approx(sfa.collisionRegion[1])
    assert cfg.helixCutTolerance == _approx(sfa.helixCutTolerance)
    assert cfg.sigmaScattering == _approx(sfa.sigmaScattering)
    assert cfg.radLengthPerSeed == _approx(sfa.radLengthPerSeed)
    assert cfg.toleranceParam == _approx(sfa.toleranceParam)
    assert cfg.deltaInvHelixDiameter == _approx(sfila.deltaInvHelixDiameter)
    assert cfg.compatSeedWeight == _approx(sfila.compatSeedWeight)
    assert cfg.impactWeightFactor == _approx(sfila.impactWeightFactor)
    assert cfg.zOriginWeightFactor == _approx(sfila.zOriginWeightFactor)
    assert cfg.maxSeedsPerSpM == sfa.maxSeedsPerSpM
    assert cfg.compatSeedLimit == sfila.compatSeedLimit
    assert cfg.seedWeightIncrement == _approx(sfila.seedWeightIncrement)
    assert cfg.numSeedIncrement == _approx(sfila.numSeedIncrement)
    assert cfg.seedConfirmation == sfa.seedConfirmation
    assert cfg.maxSeedsPerSpMConf == sfila.maxSeedsPerSpMConf
    assert cfg.maxQualitySeedsPerSpMConf == sfila.maxQualitySeedsPerSpMConf
    assert cfg.useDeltaRinsteadOfTopRadius == sfila.useDeltaRorTopRadius


def test_grid_triplet_seeding_forwards_all_applicable_fields():
    sfa = _seedFinderConfigArg(zBinsCustomLooping=[1, 2])
    sfila = _seedFilterConfigArg()
    saca = _seedingAlgorithmConfigArg()
    spga = SpacePointGridConfigArg(
        zBinEdges=[0.0, 50.0, 100.0],
        phiBinDeflectionCoverage=2,
        maxPhiBins=500,
        phi=(-1.0, 1.0),
    )
    sfoa = SeedFinderOptionsArg(bFieldInZ=1.5 * u.T)

    seq = _CapturingSequencer()
    addGridTripletSeeding(seq, "sp", saca, sfa, sfoa, sfila, spga)
    cfg = seq.algorithm.config

    _assertCommonFields(cfg, sfa, sfoa, sfila)

    # Grid/bin-finder fields, only meaningful for the grid-based algorithm.
    assert cfg.phiMin == spga.phi[0]
    assert cfg.phiMax == spga.phi[1]
    assert cfg.phiBinDeflectionCoverage == spga.phiBinDeflectionCoverage
    assert cfg.maxPhiBins == spga.maxPhiBins
    assert list(cfg.zBinEdges) == spga.zBinEdges
    assert list(cfg.zBinsCustomLooping) == sfa.zBinsCustomLooping
    assert cfg.numPhiNeighbors == saca.numPhiNeighbors
    assert list(cfg.zBinNeighborsTop) == saca.zBinNeighborsTop
    assert list(cfg.zBinNeighborsBottom) == saca.zBinNeighborsBottom
    assert cfg.useExtraCuts == saca.useExtraCuts


def test_grid_triplet_seeding_zbinedges_falls_back_to_seed_finder_arg():
    # spacePointGridConfigArg.zBinEdges unset -> falls back to
    # seedFinderConfigArg.zBinEdges (both namedtuples carry this field;
    # addStandardSeeding reads it from seedFinderConfigArg only).
    sfa = _seedFinderConfigArg(zBinEdges=[0.0, 10.0, 20.0])
    sfila = _seedFilterConfigArg()
    saca = _seedingAlgorithmConfigArg()
    spga = SpacePointGridConfigArg(phi=(-1.0, 1.0))
    sfoa = SeedFinderOptionsArg(bFieldInZ=1.5 * u.T)

    seq = _CapturingSequencer()
    addGridTripletSeeding(seq, "sp", saca, sfa, sfoa, sfila, spga)

    assert list(seq.algorithm.config.zBinEdges) == sfa.zBinEdges


def test_grid_triplet_seeding_deltazmin_defaults_to_minus_deltazmax():
    # deltaZ is a signed cut; deltaZMax alone historically meant the
    # symmetric "maximum |deltaZ|", so deltaZMin must default to
    # -deltaZMax, not the C++ default of -inf, when only deltaZMax is set.
    sfa = _seedFinderConfigArg(deltaZMin=None, deltaZMax=123 * u.mm)
    sfila = _seedFilterConfigArg()
    saca = _seedingAlgorithmConfigArg()
    spga = SpacePointGridConfigArg(phi=(-1.0, 1.0))
    sfoa = SeedFinderOptionsArg(bFieldInZ=1.5 * u.T)

    seq = _CapturingSequencer()
    addGridTripletSeeding(seq, "sp", saca, sfa, sfoa, sfila, spga)

    assert seq.algorithm.config.deltaZMin == -123 * u.mm
    assert seq.algorithm.config.deltaZMax == 123 * u.mm


def test_orthogonal_triplet_seeding_forwards_all_applicable_fields():
    sfa = _seedFinderConfigArg()
    sfila = _seedFilterConfigArg()
    saca = _seedingAlgorithmConfigArg()
    # OrthogonalTripletSeedingAlgorithm is KD-tree based: it has no
    # grid/bin-finder Config fields, so only phi is meaningful here.
    spga = SpacePointGridConfigArg(phi=(-1.0, 1.0))
    sfoa = SeedFinderOptionsArg(bFieldInZ=1.5 * u.T)

    seq = _CapturingSequencer()
    addOrthogonalTripletSeeding(seq, "sp", saca, sfa, sfoa, sfila, spga)
    cfg = seq.algorithm.config

    _assertCommonFields(cfg, sfa, sfoa, sfila)
    assert cfg.phiMin == spga.phi[0]
    assert cfg.phiMax == spga.phi[1]
    assert cfg.useExtraCuts == saca.useExtraCuts


def test_orthogonal_triplet_seeding_deltazmin_defaults_to_minus_deltazmax():
    sfa = _seedFinderConfigArg(deltaZMin=None, deltaZMax=321 * u.mm)
    sfila = _seedFilterConfigArg()
    saca = _seedingAlgorithmConfigArg()
    spga = SpacePointGridConfigArg(phi=(-1.0, 1.0))
    sfoa = SeedFinderOptionsArg(bFieldInZ=1.5 * u.T)

    seq = _CapturingSequencer()
    addOrthogonalTripletSeeding(seq, "sp", saca, sfa, sfoa, sfila, spga)

    assert seq.algorithm.config.deltaZMin == -321 * u.mm
    assert seq.algorithm.config.deltaZMax == 321 * u.mm
