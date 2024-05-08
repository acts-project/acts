#!/usr/bin/env python3

import os
import argparse
import pathlib

import acts
import acts.examples
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
    addParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    CkfConfig,
    addCKFTracks,
    TrackSelectorConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addAmbiguityResolutionML,
    AmbiguityResolutionMLConfig,
    addScoreBasedAmbiguityResolution,
    ScoreBasedAmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    addSeedFilterML,
    SeedFilterMLDBScanConfig,
)
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

u = acts.UnitConstants


parser = argparse.ArgumentParser(description="Full chain with the OpenDataDetector")
parser.add_argument(
    "--output",
    "-o",
    help="Output directory",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "odd_output",
)
parser.add_argument("--events", "-n", help="Number of events", type=int, default=100)
parser.add_argument("--skip", "-s", help="Number of events", type=int, default=0)
parser.add_argument("--edm4hep", help="Use edm4hep inputs", type=pathlib.Path)
parser.add_argument(
    "--geant4", help="Use Geant4 instead of fatras", action="store_true"
)
parser.add_argument(
    "--ttbar",
    help="Use Pythia8 (ttbar, pile-up 200) instead of particle gun",
    action="store_true",
)
parser.add_argument(
    "--ttbar-pu",
    help="Number of pile-up events for ttbar",
    type=int,
    default=200,
)
parser.add_argument(
    "--gun-multiplicity",
    help="Multiplicity of the particle gun",
    type=int,
    default=200,
)
parser.add_argument(
    "--ambi-solver",
    help="Set which ambiguity solver to use, default is the classical one",
    type=str,
    choices=["greedy", "scoring", "ML"],
    default="greedy",
)
parser.add_argument(
    "--ambi-config",
    help="Set the configuration file for the Score Based ambiguity resolution",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "ambi_config.json",
)

parser.add_argument(
    "--MLSeedFilter",
    help="Use the Ml seed filter to select seed after the seeding step",
    action="store_true",
)

args = vars(parser.parse_args())

outputDir = args["output"]
ttbar = args["ttbar"]
g4_simulation = args["geant4"]
ambi_ML = args["ambi_solver"] == "ML"
ambi_scoring = args["ambi_solver"] == "scoring"
ambi_config = args["ambi_config"]
seedFilter_ML = args["MLSeedFilter"]
geoDir = getOpenDataDetectorDirectory()
# acts.examples.dump_args_calls(locals())  # show python binding calls

oddMaterialMap = geoDir / "data/odd-material-maps.root"
oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = geoDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector, trackingGeometry, decorators = getOpenDataDetector(
    odd_dir=geoDir, mdecorator=oddMaterialDeco
)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=args["events"],
    skip=args["skip"],
    numThreads=1 if g4_simulation else -1,
    outputDir=str(outputDir),
)

if args["edm4hep"]:
    import acts.examples.edm4hep

    edm4hepReader = acts.examples.edm4hep.EDM4hepReader(
        inputPath=str(args["edm4hep"]),
        inputSimHits=[
            "PixelBarrelReadout",
            "PixelEndcapReadout",
            "ShortStripBarrelReadout",
            "ShortStripEndcapReadout",
            "LongStripBarrelReadout",
            "LongStripEndcapReadout",
        ],
        outputParticlesGenerator="particles_input",
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        graphvizOutput="graphviz",
        dd4hepDetector=detector,
        trackingGeometry=trackingGeometry,
        sortSimHitsInTime=True,
        level=acts.logging.INFO,
    )
    s.addReader(edm4hepReader)
    s.addWhiteboardAlias("particles", edm4hepReader.config.outputParticlesGenerator)

    addParticleSelection(
        s,
        config=ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            pt=(150 * u.MeV, None),
            removeNeutral=True,
        ),
        inputParticles="particles",
        outputParticles="particles_selected",
    )
else:
    if not ttbar:
        addParticleGun(
            s,
            MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
            EtaConfig(-3.0, 3.0),
            PhiConfig(0.0, 360.0 * u.degree),
            ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns
                ),
            ),
            multiplicity=args["gun_multiplicity"],
            rnd=rnd,
        )
    else:
        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=args["ttbar_pu"],
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
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
                rho=(0.0, 24 * u.mm),
                absZ=(0.0, 1.0 * u.m),
                eta=(-3.0, 3.0),
                pt=(150 * u.MeV, None),
                removeNeutral=True,
            ),
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
            rnd=rnd,
            killVolume=trackingGeometry.worldVolume,
            killAfterTime=25 * u.ns,
        )
    else:
        addFatras(
            s,
            trackingGeometry,
            field,
            preSelectParticles=(
                ParticleSelectorConfig(
                    rho=(0.0, 24 * u.mm),
                    absZ=(0.0, 1.0 * u.m),
                    eta=(-3.0, 3.0),
                    pt=(150 * u.MeV, None),
                    removeNeutral=True,
                )
                if ttbar
                else ParticleSelectorConfig()
            ),
            enableInteractions=True,
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
    (
        TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-3.0, 3.0), nHits=(9, None))
        if ttbar
        else TruthSeedRanges()
    ),
    geoSelectionConfigFile=oddSeedingSel,
    outputDirRoot=outputDir,
    # outputDirCsv=outputDir,
)

if seedFilter_ML:
    addSeedFilterML(
        s,
        SeedFilterMLDBScanConfig(
            epsilonDBScan=0.03, minPointsDBScan=2, minSeedScore=0.1
        ),
        onnxModelFile=os.path.dirname(__file__)
        + "/MLAmbiguityResolution/seedDuplicateClassifier.onnx",
        outputDirRoot=outputDir,
        # outputDirCsv=outputDir,
    )

addCKFTracks(
    s,
    trackingGeometry,
    field,
    TrackSelectorConfig(
        pt=(1.0 * u.GeV if ttbar else 0.0, None),
        absEta=(None, 3.0),
        loc0=(-4.0 * u.mm, 4.0 * u.mm),
        nMeasurementsMin=7,
        maxHoles=2,
        maxOutliers=2,
    ),
    CkfConfig(
        seedDeduplication=True,
        stayOnSeed=True,
    ),
    outputDirRoot=outputDir,
    writeCovMat=True,
    # outputDirCsv=outputDir,
)

if ambi_ML:
    addAmbiguityResolutionML(
        s,
        AmbiguityResolutionMLConfig(
            maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
        ),
        outputDirRoot=outputDir,
        # outputDirCsv=outputDir,
        onnxModelFile=os.path.dirname(__file__)
        + "/MLAmbiguityResolution/duplicateClassifier.onnx",
    )

elif ambi_scoring:
    addScoreBasedAmbiguityResolution(
        s,
        ScoreBasedAmbiguityResolutionConfig(
            minScore=0,
            minScoreSharedTracks=1,
            maxShared=2,
            maxSharedTracksPerMeasurement=2,
            pTMax=1400,
            pTMin=0.5,
            phiMax=3.14,
            phiMin=-3.14,
            etaMax=4,
            etaMin=-4,
            useAmbiguityFunction=False,
        ),
        outputDirRoot=outputDir,
        ambiVolumeFile=ambi_config,
        writeCovMat=True,
        # outputDirCsv=outputDir,
    )
else:
    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
            maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
        ),
        outputDirRoot=outputDir,
        writeCovMat=True,
        # outputDirCsv=outputDir,
    )

addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.Iterative,
    outputDirRoot=outputDir,
)

s.run()
