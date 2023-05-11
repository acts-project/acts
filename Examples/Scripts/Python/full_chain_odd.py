#!/usr/bin/env python3
import os, argparse, pathlib, contextlib, acts, acts.examples
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    addGeant4,
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
    addAmbiguityResolutionML,
    AmbiguityResolutionMLConfig,
    addVertexFitting,
    VertexFinder,
)
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector

parser = argparse.ArgumentParser(description="Full chain with the OpenDataDetector")

parser.add_argument("--events", "-n", help="Number of events", type=int, default=100)
parser.add_argument(
    "--geant4", help="Use Geant4 instead of fatras", action="store_true"
)
parser.add_argument(
    "--ttbar",
    help="Use Pythia8 (ttbar, pile-up 200) instead of particle gun",
    action="store_true",
)
parser.add_argument(
    "--MLSolver",
    help="Use the Ml Ambiguity Solver instead of the classical one",
    action="store_true",
)

args = vars(parser.parse_args())

ttbar = args["ttbar"]
g4_simulation = args["geant4"]
ambiguity_MLSolver = args["MLSolver"]
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

# TODO Geant4 currently crashes with FPE monitoring
with acts.FpeMonitor() if not g4_simulation else contextlib.nullcontext():
    s = acts.examples.Sequencer(
        events=args["events"],
        numThreads=1,
        outputDir=str(outputDir),
    )

    if not ttbar:
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
            npileup=50,
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            rnd=rnd,
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
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
            preSelectParticles=ParticleSelectorConfig(
                eta=(-3.0, 3.0),
                absZ=(0, 1e4),
                rho=(0, 1e3),
                pt=(150 * u.MeV, None),
                removeNeutral=True,
            ),
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
            rnd=rnd,
            killVolume=acts.Volume.makeCylinderVolume(r=1100, halfZ=3000),
        )
    else:
        addFatras(
            s,
            trackingGeometry,
            field,
            preSelectParticles=ParticleSelectorConfig(
                eta=(-3.0, 3.0),
                pt=(150 * u.MeV, None),
                removeNeutral=True,
            )
            if ttbar
            else ParticleSelectorConfig(),
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
            rnd=rnd,
        )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=oddDigiConfig,
        outputDirRoot=outputDir,
        # outputDirCsv=outputDir,
        rnd=rnd,
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-3.0, 3.0), nHits=(9, None))
        if ttbar
        else TruthSeedRanges(),
        geoSelectionConfigFile=oddSeedingSel,
        outputDirRoot=outputDir,
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFPerformanceConfig(
            ptMin=1.0 * u.GeV if ttbar else 0.0,
            nMeasurementsMin=7,
        ),
        TrackSelectorRanges(
            pt=(1.0 * u.GeV, None),
            absEta=(None, 3.0),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
        ),
        outputDirRoot=outputDir,
        # outputDirCsv=outputDir,
    )

    if ambiguity_MLSolver:
        addAmbiguityResolutionML(
            s,
            AmbiguityResolutionMLConfig(nMeasurementsMin=7),
            CKFPerformanceConfig(
                ptMin=1.0 * u.GeV if ttbar else 0.0, nMeasurementsMin=7
            ),
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
            onnxModelFile=os.path.dirname(__file__)
            + "/MLAmbiguityResolution/duplicateClassifier.onnx",
        )
    else:
        addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(
                maximumSharedHits=3, maximumIterations=10000, nMeasurementsMin=7
            ),
            CKFPerformanceConfig(
                ptMin=1.0 * u.GeV if ttbar else 0.0,
                nMeasurementsMin=7,
            ),
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
        )

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=outputDir,
    )

    s.run()
