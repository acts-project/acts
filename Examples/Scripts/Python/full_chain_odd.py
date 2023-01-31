#!/usr/bin/env python3
import pathlib, acts, acts.examples
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    addCKFTracks,
    CKFPerformanceConfig,
    TrackSelectorRanges,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
)
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector

ttbar_pu200 = False
g4_simulation = False
u = acts.UnitConstants
geoDir = getOpenDataDetectorDirectory()
outputDir = pathlib.Path.cwd() / "odd_output"
# acts.examples.dump_args_calls(locals())  # show python binding calls

oddMaterialMap = geoDir / "data/odd-material-maps.root"
oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = geoDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector, trackingGeometry, decorators = getOpenDataDetector(
    geoDir, mdecorator=oddMaterialDeco
)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

with acts.FpeMonitor():
    s = acts.examples.Sequencer(events=100, numThreads=1, outputDir=str(outputDir))

    if not ttbar_pu200:
        addParticleGun(
            s,
            MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
            EtaConfig(-3.0, 3.0, uniform=True),
            ParticleConfig(2, acts.PdgParticle.eMuon, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=50,
            rnd=rnd,
        )
    else:
        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=200,
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            rnd=rnd,
            outputDirRoot=outputDir,
        )
    if g4_simulation:
        if s.config.numThreads != 1:
            raise ValueError("Geant 4 simulation does not support multi-threading")

        # Pythia can sometime simulate particles outside the world volume, a cut on the Z of the track help mitigate this effect
        # Older version of G4 might not work, this as has been tested on version `geant4-11-00-patch-03`
        # For more detail see issue #1578
        addGeant4(
            s,
            detector,
            trackingGeometry,
            field,
            preselectParticles=ParticleSelectorConfig(
                eta=(-3.0, 3.0),
                absZ=(0, 1e4),
                pt=(150 * u.MeV, None),
                removeNeutral=True,
            ),
            outputDirCsv=outputDir,
            rnd=rnd,
        )
    else:
        addFatras(
            s,
            trackingGeometry,
            field,
            ParticleSelectorConfig(
                eta=(-3.0, 3.0), pt=(150 * u.MeV, None), removeNeutral=True
            )
            if ttbar_pu200
            else ParticleSelectorConfig(),
            outputDirRoot=outputDir,
            rnd=rnd,
        )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=oddDigiConfig,
        outputDirRoot=outputDir,
        rnd=rnd,
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-3.0, 3.0), nHits=(9, None))
        if ttbar_pu200
        else TruthSeedRanges(),
        geoSelectionConfigFile=oddSeedingSel,
        outputDirRoot=outputDir,
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFPerformanceConfig(
            ptMin=1.0 * u.GeV if ttbar_pu200 else 0.0, nMeasurementsMin=6
        ),
        TrackSelectorRanges(
            pt=(1.0 * u.GeV, None),
            absEta=(None, 3.0),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
            removeNeutral=True,
        ),
        outputDirRoot=outputDir,
    )

    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(maximumSharedHits=3),
        CKFPerformanceConfig(
            ptMin=1.0 * u.GeV if ttbar_pu200 else 0.0, nMeasurementsMin=6
        ),
        outputDirRoot=outputDir,
    )

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=outputDir,
    )

    s.run()
