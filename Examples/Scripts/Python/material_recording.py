#!/usr/bin/env python3
import os
import warnings
from pathlib import Path

import acts
from acts.examples import (
    GaussianVertexGenerator,
    ParametricParticleGenerator,
    FixedMultiplicityGenerator,
    EventGenerator,
    RandomNumbers,
)

import acts.examples.dd4hep
import acts.examples.geant4
import acts.examples.geant4.dd4hep
from common import getOpenDataDetectorDirectory

u = acts.UnitConstants

_material_recording_executed = False


def runMaterialRecording(g4geo, outputDir, tracksPerEvent=10000, s=None):
    global _material_recording_executed
    if _material_recording_executed:
        warnings.warn("Material recording already ran in this process. Expect crashes")
    _material_recording_executed = True

    rnd = RandomNumbers(seed=228)

    s = s or acts.examples.Sequencer(events=1000, numThreads=1)

    evGen = EventGenerator(
        level=acts.logging.INFO,
        generators=[
            EventGenerator.Generator(
                multiplicity=FixedMultiplicityGenerator(n=1),
                vertex=GaussianVertexGenerator(
                    stddev=acts.Vector4(0, 0, 0, 0),
                    mean=acts.Vector4(0, 0, 0, 0),
                ),
                particles=ParametricParticleGenerator(
                    p=(1 * u.GeV, 10 * u.GeV),
                    eta=(-4, 4),
                    numParticles=tracksPerEvent,
                    etaUniform=True,
                ),
            )
        ],
        outputParticles="particles_initial",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    g4AlgCfg = acts.examples.geant4.materialRecordingConfig(
        level=acts.logging.INFO,
        detector=g4geo,
        inputParticles=evGen.config.outputParticles,
        outputMaterialTracks="material_tracks",
        randomNumbers=rnd,
    )

    g4AlgCfg.detectorConstruction = g4geo

    g4Alg = acts.examples.geant4.Geant4Simulation(
        level=acts.logging.INFO, config=g4AlgCfg
    )

    s.addAlgorithm(g4Alg)

    s.addWriter(
        acts.examples.RootMaterialTrackWriter(
            prePostStep=True,
            recalculateTotals=True,
            collection="material_tracks",
            filePath=os.path.join(outputDir, "geant4_material_tracks.root"),
            level=acts.logging.INFO,
        )
    )

    return s


if "__main__" == __name__:

    dd4hepSvc = acts.examples.dd4hep.DD4hepGeometryService(
        xmlFileNames=[str(getOpenDataDetectorDirectory() / "xml/OpenDataDetector.xml")]
    )
    g4geo = acts.examples.geant4.dd4hep.DDG4DetectorConstruction(dd4hepSvc)

    runMaterialRecording(g4geo=g4geo, tracksPerEvent=100, outputDir=os.getcwd()).run()
