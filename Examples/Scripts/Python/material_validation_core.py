#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

import acts
from acts import (
    MaterialValidater,
    IntersectionMaterialAssigner,
    MaterialMapJsonConverter,
)

from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    RootMaterialTrackWriter,
    MaterialValidation,
)

from acts.examples.dd4hep import (
    DD4hepDetector,
    DD4hepDetectorOptions,
    DD4hepGeometryService,
)

from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory


def runMaterialValidation(s, ntracks, surfaces, outputFile, seed, loglevel):
    # IO for material tracks reading
    wb = WhiteBoard(acts.logging.INFO)

    rnd = acts.examples.RandomNumbers(seed=seed)

    # Assignment setup : Intersection assigner
    materialAssingerConfig = IntersectionMaterialAssigner.Config()
    materialAssingerConfig.surfaces = surfaces
    materialAssinger = IntersectionMaterialAssigner(materialAssingerConfig, loglevel)

    # Validater setup
    materialValidaterConfig = MaterialValidater.Config()
    materialValidaterConfig.materialAssigner = materialAssinger
    materialValidater = MaterialValidater(materialValidaterConfig, loglevel)

    # Validation Algorithm
    materialValidationConfig = MaterialValidation.Config()
    materialValidationConfig.materialValidater = materialValidater
    materialValidationConfig.outputMaterialTracks = "recorded-material-tracks"
    materialValidationConfig.ntracks = ntracks
    materialValidationConfig.randomNumberSvc = rnd
    materialValidation = MaterialValidation(materialValidationConfig, loglevel)
    s.addAlgorithm(materialValidation)

    # Add the mapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=materialValidationConfig.outputMaterialTracks,
            filePath=outputFile + ".root",
            storeSurface=True,
            storeVolume=True,
        )
    )
    # return the sequencer
    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser()

    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to process"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=100, help="Number of tracks per event"
    )
    p.add_argument(
        "-j", "--threads", type=int, default=-1, help="Number of threads in parallel"
    )
    p.add_argument(
        "-m", "--map", type=str, default="", help="Input file for the material map"
    )
    p.add_argument("-o", "--output", type=str, default="", help="Output file name")

    p.add_argument(
        "--experimental",
        action=argparse.BooleanOptionalAction,
        help="Construct experimental geometry",
    )

    args = p.parse_args()

    decorators = None
    if args.map != "":
        decorators = acts.IMaterialDecorator.fromFile(args.map)

    if args.experimental:
        odd_xml = getOpenDataDetectorDirectory() / "xml" / "OpenDataDetector.xml"

        # Create the dd4hep geometry service and detector
        dd4hepConfig = DD4hepGeometryService.Config()
        dd4hepConfig.logLevel = acts.logging.INFO
        dd4hepConfig.xmlFileNames = [str(odd_xml)]
        dd4hepGeometryService = DD4hepGeometryService(dd4hepConfig)
        dd4hepDetector = DD4hepDetector(dd4hepGeometryService)

        cOptions = DD4hepDetectorOptions(logLevel=acts.logging.INFO, emulateToGraph="")
        cOptions.materialDecorator = decorators

        # Context and options
        geoContext = acts.GeometryContext()
        [detector, contextors, store] = dd4hepDetector.finalize(geoContext, cOptions)

        materialSurfaces = detector.extractMaterialSurfaces()

    else:
        [detector, trackingGeometry, decorators] = getOpenDataDetector(decorators)

        materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    s = acts.examples.Sequencer(events=args.events, numThreads=args.threads)

    runMaterialValidation(
        s, args.tracks, materialSurfaces, args.output, 42, acts.logging.INFO
    ).run()
