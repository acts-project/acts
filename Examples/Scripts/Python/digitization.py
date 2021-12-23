#!/usr/bin/env python3
import os
from pathlib import Path
from typing import Optional

import acts
import acts.examples


u = acts.UnitConstants


def configureDigitization(
    trackingGeometry,
    field,
    outputDir: Path,
    particlesInput: Optional[Path] = None,
    outputRoot=True,
    outputCsv=True,
    s=None,
):

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    csv_dir = os.path.join(outputDir, "csv")
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)

    # Input
    rnd = acts.examples.RandomNumbers(seed=42)

    particleCollection = "particles_input"
    if particlesInput is None:
        evGen = acts.examples.EventGenerator(
            level=acts.logging.INFO,
            generators=[
                acts.examples.EventGenerator.Generator(
                    multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
                    vertex=acts.examples.GaussianVertexGenerator(
                        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                    ),
                    particles=acts.examples.ParametricParticleGenerator(
                        p=(1 * u.GeV, 10 * u.GeV),
                        eta=(-2, 2),
                        phi=(0, 360 * u.degree),
                        randomizeCharge=True,
                        numParticles=4,
                    ),
                )
            ],
            outputParticles=particleCollection,
            randomNumbers=rnd,
        )
    else:
        # Read input from input collection (e.g. Pythia8 output)
        evGen = acts.examples.RootParticleReader(
            level=acts.logging.INFO,
            particleCollection=particleCollection,
            filePath=str(particlesInput),
            orderedEvents=False,
        )

    # Simulation
    simAlg = acts.examples.FatrasSimulation(
        level=acts.logging.INFO,
        inputParticles=particleCollection,
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
    )

    # Digitization
    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(
            str(
                srcdir
                / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
            )
        ),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits=simAlg.config.outputSimHits,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )

    s.addReader(evGen)
    s.addAlgorithm(simAlg)
    s.addAlgorithm(digiAlg)

    if outputRoot:
        rmwConfig = acts.examples.RootMeasurementWriter.Config(
            inputMeasurements=digiAlg.config.outputMeasurements,
            inputClusters=digiAlg.config.outputClusters,
            inputSimHits=simAlg.config.outputSimHits,
            inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
            filePath=str(outputDir / f"{digiAlg.config.outputMeasurements}.root"),
            trackingGeometry=trackingGeometry,
        )
        rmwConfig.addBoundIndicesFromDigiConfig(digiAlg.config)
        s.addWriter(acts.examples.RootMeasurementWriter(rmwConfig, acts.logging.INFO))

    if outputCsv:
        csv_dir = outputDir / "csv"
        csv_dir.mkdir(parents=True, exist_ok=True)
        s.addWriter(
            acts.examples.CsvMeasurementWriter(
                level=acts.logging.VERBOSE,
                inputMeasurements=digiAlg.config.outputMeasurements,
                inputClusters=digiAlg.config.outputClusters,
                inputSimHits=simAlg.config.outputSimHits,
                inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
                outputDir=str(csv_dir),
            )
        )

    return s


if "__main__" == __name__:
    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    configureDigitization(trackingGeometry, field, outputDir=Path.cwd()).run()
