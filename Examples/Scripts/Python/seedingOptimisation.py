#!/usr/bin/env python3
import os
import dataclasses
import functools

import acts
import acts.examples

from common import addPythia8

u = acts.UnitConstants

def runSimulation(trackingGeometry, field, rnd, outputDir):

    csv_dir = os.path.join(outputDir, "csv")
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)

# Sequencer
    s = acts.examples.Sequencer(
        events=1000, numThreads=-1, logLevel=acts.logging.INFO
    )

    evGen = addPythia8(s, rnd, hardProcess = ["Top:qqbar2ttbar=on"])

    # Simulation
    simAlg = acts.examples.FatrasSimulation(
        level=acts.logging.INFO,
        inputParticles=evGen.config.outputParticles,
        # inputParticles=evGen.config.particleCollection,
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
    )

    s.addAlgorithm(simAlg)

    # Output
    s.addWriter(
        acts.examples.CsvParticleWriter(
            level=acts.logging.INFO,
            outputDir=outputDir + "/csv",
            inputParticles="particles_final",
            outputStem="particles_final",
        )
    )

    s.addWriter(
        acts.examples.CsvSimHitWriter(
            level=acts.logging.INFO,
            inputSimHits="simhits",
            outputDir=outputDir + "/csv",
            outputStem="hits",
        )
    )

    s.run()

    return()


def runSeeding(trackingGeometry, field, rnd, outputDir,  gridConfig, seedFilterConfig, seedFinderConfig):


    partReader = acts.examples.CsvParticleReader(
        level=acts.logging.INFO,
        inputDir=outputDir + "/csv",
        inputStem="particles_final",
        outputParticles="particles_final",
    )

    hitReader = acts.examples.CsvSimHitReader(
        level=acts.logging.INFO,
        inputDir=outputDir + "/csv",
        inputStem="hits",
        outputSimHits="simhits",
    )

    # Digitization
    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(
            "../Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits=hitReader.config.outputSimHits,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)

    selAlg = acts.examples.TruthSeedSelector(
        level=acts.logging.INFO,
        ptMin=1 * u.GeV,
        eta=(-2.5, 2.5),
        nHitsMin=9,
        inputParticles=partReader.config.outputParticles,
        inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
        outputParticles="particles_selected",
    )

    inputParticles = selAlg.config.outputParticles

    spAlg = acts.examples.SpacePointMaker(
        level=acts.logging.INFO,
        inputSourceLinks=digiCfg.outputSourceLinks,
        inputMeasurements=digiCfg.outputMeasurements,
        outputSpacePoints="spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.readJsonGeometryList(
            "../Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=acts.logging.INFO,
        inputSpacePoints=[spAlg.config.outputSpacePoints],
        outputSeeds="seeds",
        outputProtoTracks="prototracks",
        gridConfig=gridConfig,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
    )

    parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
        level=acts.logging.INFO,
        inputProtoTracks=seedingAlg.config.outputProtoTracks,
        inputSpacePoints=[spAlg.config.outputSpacePoints],
        inputSourceLinks=digiCfg.outputSourceLinks,
        outputTrackParameters="estimatedparameters",
        outputProtoTracks="prototracks_estimated",
        trackingGeometry=trackingGeometry,
        magneticField=field,
    )

    s = acts.examples.Sequencer(
        events=1000,
        numThreads=-1,
        logLevel=acts.logging.INFO,
    )

    s.addReader(partReader)
    s.addReader(hitReader)
    s.addAlgorithm(digiAlg)
    s.addAlgorithm(selAlg)
    s.addAlgorithm(spAlg)
    s.addAlgorithm(seedingAlg)
    s.addAlgorithm(parEstimateAlg)

    seedingPerformaces = acts.examples.SeedingPerformanceWriter(
        level=acts.logging.INFO,
        inputProtoTracks=seedingAlg.config.outputProtoTracks,
        inputParticles=inputParticles,
        inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
        filePath=outputDir + "/performance_seeding_hists.root",
    )

    s.addWriter(seedingPerformaces)
    s.run()

    if seedingPerformaces.totalParticles != 0:
        efficiency = seedingPerformaces.totalMatchedParticles / seedingPerformaces.totalParticles
    else:
        efficiency = 0
        
    if seedingPerformaces.totalSeeds != 0:
        fakeRate = (seedingPerformaces.totalSeeds - seedingPerformaces.totalMatchedSeeds)  / seedingPerformaces.totalSeeds
    else:
        fakeRate = 1

    if seedingPerformaces.totalMatchedParticles != 0:
        duplicateRate = seedingPerformaces.totalDuplicatedParticles / seedingPerformaces.totalMatchedParticles
    else:
        duplicateRate = 1
        
    if seedingPerformaces.totalMatchedParticles != 0:
        avgDuplicate = (seedingPerformaces.totalMatchedSeeds - seedingPerformaces.totalMatchedParticles)  / seedingPerformaces.totalMatchedParticles
    else:
        avgDuplicate = 0        

    return (
        seedingPerformaces.totalSeeds,
        seedingPerformaces.totalMatchedSeeds,
        seedingPerformaces.totalParticles,
        seedingPerformaces.totalMatchedParticles,
        seedingPerformaces.totalDuplicatedParticles,
        efficiency,
        fakeRate,
        duplicateRate,
        avgDuplicate,
    )


class Attribute:
    def __init__(self, _type, min, max):
        assert isinstance(min, _type) and isinstance(max, _type)
        assert type(min) == type(max), f"{type(min)}, {type(max)}"
        self.min = min
        self.max = max


limits = {"maxSeedsPerSpM": Attribute(int, 1, 10)}


# @dataclasses.dataclass
# class SeedingConfig:
#     maxSeedsPerSpM: int = 1
#     rMax: float = 200 * u.mm
#     deltaRMin: float = 1 * u.mm
#     deltaRMax: float = 60 * u.mm
#     sigmaScattering: float = 10
#     impactMax: float = 3 * u.mm
#     minPt: float = 500 * u.MeV

if "__main__" == __name__:
    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    rnd = acts.examples.RandomNumbers(seed=42)

    outputDir="seeding"
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)


    runSimulation(trackingGeometry, field, rnd, outputDir)

    from orion.client import build_experiment

    storage = {
        "database": {
            "type": "ephemeraldb",
        },
    }

    space = {
        "maxSeedsPerSpM": "uniform(1, 10)",
        "minPt": "uniform(100, 1000)",
        "deltaRMax": "uniform(10, 100)",
        "deltaRMin": "uniform(1, 10)",
        "radLengthPerSeed": "uniform(0.01, 0.1)",
        "compatSeedWeight": "uniform(100, 1000)",
        "impactWeightFactor": "uniform(0.5, 5)",
    }

    experiment = build_experiment(
        "random-exp",
        space=space,
        storage=storage,
        # algorithms={"tpe": {"n_initial_points": 5}},
    )

    def evaluate(maxSeedsPerSpM, minPt, deltaRMax, deltaRMin, radLengthPerSeed, compatSeedWeight, impactWeightFactor):
        print("evaluate")

        if deltaRMax < deltaRMin:
            deltaRMin, deltaRMax = deltaRMax, deltaRMin

        gridConfig = acts.SpacePointGridConfig(
            bFieldInZ=1.99724 * u.T,
            minPt=minPt * u.MeV,
            rMax=200 * u.mm,
            zMax=2000 * u.mm,
            zMin=-2000 * u.mm,
            deltaRMax=deltaRMax * u.mm,
            # cotThetaMax=7.40627,  # 2.7 eta
            cotThetaMax=27.29,  # 4 eta
        )

        seedFilterConfig = acts.SeedFilterConfig(
            maxSeedsPerSpM=round(maxSeedsPerSpM), deltaRMin=deltaRMin * u.mm, compatSeedWeight=compatSeedWeight, impactWeightFactor=impactWeightFactor
        )
        seedFinderConfig = acts.SeedfinderConfig(
            rMax=gridConfig.rMax,
            deltaRMin=seedFilterConfig.deltaRMin,
            deltaRMax=gridConfig.deltaRMax,
            collisionRegionMin=-250 * u.mm,
            collisionRegionMax=250 * u.mm,
            zMin=gridConfig.zMin,
            zMax=gridConfig.zMax,
            cotThetaMax=gridConfig.cotThetaMax,
            sigmaScattering=5,
            radLengthPerSeed=radLengthPerSeed,
            minPt=gridConfig.minPt,
            bFieldInZ=gridConfig.bFieldInZ,
            beamPos=acts.Vector2(0 * u.mm, 0 * u.mm),
            impactMax=3*u.mm,
        )

        (
            nTotalSeeds,
            nTotalMatchedSeeds,
            nTotalParticles,
            nTotalMatchedParticles,
            nTotalDuplicatedParticles,
            efficiency,
            fakeRate,
            duplicateRate,
            avgDuplicate,
        ) = runSeeding(trackingGeometry, field, rnd, outputDir,  gridConfig, seedFilterConfig, seedFinderConfig)

        K = 1000
        effScore = efficiency - (100 * fakeRate + avgDuplicate) / K
        print("efficiency : ", efficiency, "fakeRate : ", fakeRate, "duplicateRate : ",duplicateRate, "effScore : ",effScore)

        objective = 1 - effScore

        return [
            {"name": "score", "type": "objective", "value": objective},
        ]


    print("begin workon")
    experiment.workon(evaluate, max_trials=50)
    print("workon done")

    experiment.plot.regret().show()

    df = experiment.to_pandas()

    best = df.iloc[df.objective.idxmin()]
    print(best)
