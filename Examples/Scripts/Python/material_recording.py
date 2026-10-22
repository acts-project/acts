#!/usr/bin/env python3

import argparse
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
import acts.examples.hepmc3
from acts.examples.odd import getOpenDataDetector
from acts.examples.root import RootMaterialTrackWriter

u = acts.UnitConstants


def runMaterialRecording(
    detector,
    s,
    tracksPerEvent=10000,
    etaRange=(-4.0, 4.0),
    phiRange=(0.0, 360.0 * u.degree),
    materialTrackCollectionName="material_tracks",
    outputFileBase="geant4_material_tracks",
):

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
                    phi=phiRange,
                    numParticles=tracksPerEvent,
                    etaUniform=True,
                ),
            )
        ],
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    hepmc3Converter = acts.examples.hepmc3.HepMC3InputConverter(
        level=acts.logging.INFO,
        inputEvent=evGen.config.outputEvent,
        outputParticles="particles_initial",
        outputVertices="vertices_initial",
        mergePrimaries=False,
    )
    s.addAlgorithm(hepmc3Converter)

    g4Alg = acts.examples.geant4.Geant4MaterialRecording(
        level=acts.logging.INFO,
        detector=detector,
        randomNumbers=rnd,
        inputParticles=hepmc3Converter.config.outputParticles,
        outputMaterialTracks=materialTrackCollectionName,
        recordElementFractions=False,
    )

    s.addAlgorithm(g4Alg)

    s.addWriter(
        RootMaterialTrackWriter(
            prePostStep=True,
            recalculateTotals=True,
            inputMaterialTracks=materialTrackCollectionName,
            treeName=materialTrackCollectionName,
            filePath=str(outputFileBase) + ".root",
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
        "-i", "--input", type=str, default="", help="input (GDML/SQL) file (optional)"
    )
    p.add_argument(
        "--eta-range",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        default=(-4.0, 4.0),
        help="Eta range for generated particles",
    )
    p.add_argument(
        "--phi-range",
        nargs=2,
        type=float,
        metavar=("MIN_DEG", "MAX_DEG"),
        default=(0.0, 360.0),
        help="Phi range in degree for generated particles",
    )
    p.add_argument(
        "--material-track-collection",
        type=str,
        default="material_tracks",
        help="Output material track collection name",
    )
    p.add_argument(
        "-o",
        "--output",
        type=str,
        default="material_geant4",
        help="Output file stem (without extension)",
    )

    args = p.parse_args()

    detector = None
    if args.input == "":
        detector = getOpenDataDetector()
    elif args.input.endswith(".gdml"):
        detector = acts.examples.geant4.GdmlDetector(path=args.input)
    elif args.input.endswith(".sqlite") or args.input.endswith(".db"):
        gmdConfig = acts.geomodel.GeoModelDetector.Config(path=args.input)
        detector = acts.geomodel.GeoModelDetector(gmdConfig)

    runMaterialRecording(
        detector=detector,
        s=acts.examples.Sequencer(events=args.events, numThreads=1),
        tracksPerEvent=args.tracks,
        etaRange=tuple(args.eta_range),
        phiRange=(args.phi_range[0] * u.degree, args.phi_range[1] * u.degree),
        materialTrackCollectionName=args.material_track_collection,
        outputFileBase=args.output,
    ).run()


if "__main__" == __name__:
    main()
