#!/usr/bin/env python3

import os
import warnings
import argparse

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
from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants

_material_recording_executed = False


def runMaterialRecording(
    detectorConstructionFactory,
    outputDir,
    tracksPerEvent=10000,
    s=None,
    etaRange=(-4, 4),
):
    global _material_recording_executed
    if _material_recording_executed:
        warnings.warn("Material recording already ran in this process. Expect crashes")
    _material_recording_executed = True

    rnd = RandomNumbers(seed=228)

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
                    pdg=acts.PdgParticle.eInvalid,
                    charge=0,
                    randomizeCharge=False,
                    mass=0,
                    p=(1 * u.GeV, 10 * u.GeV),
                    eta=etaRange,
                    numParticles=tracksPerEvent,
                    etaUniform=True,
                ),
            )
        ],
        outputParticles="particles_initial",
        outputVertices="vertices_initial",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    g4Alg = acts.examples.geant4.Geant4MaterialRecording(
        level=acts.logging.INFO,
        detectorConstructionFactory=detectorConstructionFactory,
        randomNumbers=rnd,
        inputParticles=evGen.config.outputParticles,
        outputMaterialTracks="material_tracks",
    )

    s.addAlgorithm(g4Alg)

    s.addWriter(
        acts.examples.RootMaterialTrackWriter(
            prePostStep=True,
            recalculateTotals=True,
            inputMaterialTracks="material_tracks",
            filePath=os.path.join(outputDir, "geant4_material_tracks.root"),
            level=acts.logging.INFO,
        )
    )

    return s


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to generate"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=100, help="Particle tracks per event"
    )
    p.add_argument(
        "-i", "--input", type=str, default="", help="GDML input file (optional)"
    )

    args = p.parse_args()

    detectorConstructionFactory = None
    if args.input != "":
        detectorConstructionFactory = (
            acts.examples.geant4.GdmlDetectorConstructionFactory(args.input)
        )
    else:
        detector, trackingGeometry, decorators = getOpenDataDetector()

        detectorConstructionFactory = (
            acts.examples.geant4.dd4hep.DDG4DetectorConstructionFactory(detector)
        )

    runMaterialRecording(
        detectorConstructionFactory=detectorConstructionFactory,
        tracksPerEvent=args.tracks,
        outputDir=os.getcwd(),
        s=acts.examples.Sequencer(events=args.events, numThreads=1),
    ).run()


if "__main__" == __name__:
    main()
