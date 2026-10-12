#!/usr/bin/env python3
"""Record Geant4 material tracks through the external nODD pixel detector.

The generation and recording pipeline follows
Examples/Scripts/Python/material_recording.py, with the nODD detector adapter.
"""

import argparse
import math
from pathlib import Path

import acts
import acts.examples
import acts.examples.geant4
import acts.examples.hepmc3
from acts.examples import (
    EventGenerator,
    FixedMultiplicityGenerator,
    GaussianVertexGenerator,
    ParametricParticleGenerator,
    RandomNumbers,
)
from acts.examples.nodd import getNoddDetector
from acts.examples.root import RootMaterialTrackWriter

u = acts.UnitConstants


def runMaterialRecording(
    detector,
    s,
    tracksPerEvent=10000,
    etaRange=(-4.0, 4.0),
    phiRange=(0.0, 360.0 * u.degree),
    materialTrackCollectionName="material_tracks",
    outputFileBase="nodd_geant4_material_tracks",
    seed=228,
):
    """Add the standard neutral-geantino material recording chain to a sequencer."""
    rnd = RandomNumbers(seed=seed)
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
    converter = acts.examples.hepmc3.HepMC3InputConverter(
        level=acts.logging.INFO,
        inputEvent=evGen.config.outputEvent,
        outputParticles="particles_initial",
        outputVertices="vertices_initial",
        mergePrimaries=False,
    )
    s.addAlgorithm(converter)
    s.addAlgorithm(
        acts.examples.geant4.Geant4MaterialRecording(
            level=acts.logging.INFO,
            detector=detector,
            randomNumbers=rnd,
            inputParticles=converter.config.outputParticles,
            outputMaterialTracks=materialTrackCollectionName,
            recordElementFractions=False,
        )
    )
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


def _positive_int(value):
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodd-dir", type=Path, help="nODD checkout (or NODD_PATH)")
    parser.add_argument("--build-dir", type=Path, help="nODD build directory override")
    parser.add_argument("--compact", type=Path, help="nODD compact XML override")
    parser.add_argument("-n", "--events", type=_positive_int, default=1000)
    parser.add_argument("-t", "--tracks", type=_positive_int, default=100)
    parser.add_argument("--seed", type=int, default=228)
    parser.add_argument(
        "--eta-range", nargs=2, type=float, metavar=("MIN", "MAX"), default=(-4.0, 4.0)
    )
    parser.add_argument(
        "--phi-range",
        nargs=2,
        type=float,
        metavar=("MIN_DEG", "MAX_DEG"),
        default=(0.0, 360.0),
    )
    parser.add_argument("--material-track-collection", default="material_tracks")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("nodd_material_geant4"),
        help="Output ROOT file stem (without extension)",
    )
    args = parser.parse_args()
    for name, bounds in (("eta", args.eta_range), ("phi", args.phi_range)):
        if not all(math.isfinite(v) for v in bounds) or bounds[0] >= bounds[1]:
            parser.error(f"--{name}-range requires two finite increasing values")
    if args.seed < 0:
        parser.error("--seed must be nonnegative")

    detector = getNoddDetector(
        nodd_dir=args.nodd_dir,
        build_dir=args.build_dir,
        compact_file=args.compact,
        gen3=False,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    runMaterialRecording(
        detector=detector,
        s=acts.examples.Sequencer(events=args.events, numThreads=1),
        tracksPerEvent=args.tracks,
        etaRange=tuple(args.eta_range),
        phiRange=tuple(value * u.degree for value in args.phi_range),
        materialTrackCollectionName=args.material_track_collection,
        outputFileBase=args.output,
        seed=args.seed,
    ).run()


if __name__ == "__main__":
    main()
