#!/usr/bin/env python3

import argparse
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory

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
from acts.examples.geant4_multiprocessing import (
    mergeGeant4RootFiles,
    runGeant4Processes,
)

u = acts.UnitConstants


def runMaterialRecording(
    detector,
    s,
    tracksPerEvent=10000,
    etaRange=(-4.0, 4.0),
    phiRange=(0.0, 360.0 * u.degree),
    materialTrackCollectionName="material_tracks",
    outputFileBase="geant4_material_tracks",
    seed=228,
):

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


def _recordMaterialEventRange(
    begin, end, *, detectorFactory, outputDir, recordingOptions
):
    output = Path(outputDir) / f"events_{begin}_{end}"
    output.mkdir()
    detector = detectorFactory()
    sequencer = acts.examples.Sequencer(
        events=end - begin, skip=begin, numThreads=1, outputDir=str(output)
    )
    runMaterialRecording(
        detector=detector,
        s=sequencer,
        outputFileBase=str(output / "material"),
        **recordingOptions,
    ).run()


def runMaterialRecordingParallel(
    detectorFactory,
    *,
    events,
    processes,
    skip=0,
    outputFileBase="geant4_material_tracks",
    **recordingOptions,
):
    """Build a detector per process and merge records in global event order.

    The factory must be importable/picklable and construct the detector locally.
    Call under a main guard. Other options are forwarded to runMaterialRecording;
    every worker uses the same seed (228 by default).
    """
    output = Path(str(outputFileBase) + ".root")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="material-workers-", dir=output.parent) as temporary:
        ranges = runGeant4Processes(
            partial(
                _recordMaterialEventRange,
                detectorFactory=detectorFactory,
                outputDir=temporary,
                recordingOptions=recordingOptions,
            ),
            events=events,
            processes=processes,
            skip=skip,
        )
        mergeGeant4RootFiles(
            [
                Path(temporary) / f"events_{begin}_{end}" / "material.root"
                for begin, end in ranges
            ],
            output,
        )


def _makeDetector(path):
    if not path:
        return getOpenDataDetector()
    if path.endswith(".gdml"):
        return acts.examples.geant4.GdmlDetector(path=path)
    if path.endswith((".sqlite", ".db")):
        from acts import geomodel

        return geomodel.GeoModelDetector(geomodel.GeoModelDetector.Config(path=path))
    raise ValueError(f"Unsupported detector input: {path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to generate"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=100, help="Particle tracks per event"
    )
    p.add_argument(
        "-j", "--processes", type=int, default=1, help="Independent Geant4 processes"
    )
    p.add_argument("--skip", type=int, default=0, help="First global ACTS event number")
    p.add_argument("--seed", type=int, default=228)
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

    if args.events < 0 or args.skip < 0 or args.processes < 1:
        p.error("events and skip must be nonnegative; processes must be positive")
    detectorFactory = partial(_makeDetector, args.input)
    options = dict(
        tracksPerEvent=args.tracks,
        etaRange=tuple(args.eta_range),
        phiRange=(args.phi_range[0] * u.degree, args.phi_range[1] * u.degree),
        materialTrackCollectionName=args.material_track_collection,
        outputFileBase=args.output,
        seed=args.seed,
    )
    if args.processes == 1:
        runMaterialRecording(
            detector=detectorFactory(),
            s=acts.examples.Sequencer(events=args.events, skip=args.skip, numThreads=1),
            **options,
        ).run()
    else:
        runMaterialRecordingParallel(
            detectorFactory,
            events=args.events,
            processes=args.processes,
            skip=args.skip,
            **options,
        )


if "__main__" == __name__:
    main()
