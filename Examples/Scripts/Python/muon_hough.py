#!/usr/bin/env python3
import os
import argparse

import acts

from acts.examples import (
    CsvMuonSpacePointReader,
    CsvMuonSegmentReader,
    MuonHoughSeeder,
    RootMuonSpacePointReader,
)


u = acts.UnitConstants
rnd = acts.examples.RandomNumbers(seed=42)


def runHoughFromRoot(inFile: str, nEvents: int):
    s = acts.examples.Sequencer(
        events=nEvents, numThreads=1, logLevel=acts.logging.VERBOSE
    )

    # Read input space points from input csv files
    evReader = RootMuonSpacePointReader(
        filePath=inFile,
        outputSpacePoints="MuonSpacePoints",
        level=acts.logging.VERBOSE,
    )
    s.addReader(evReader)

    s.run()


def runHoughFromCsv(inDir: str, nEvents: int):
    # create temporary file with pixel SPs and run the seeding

    s = acts.examples.Sequencer(
        events=nEvents, numThreads=1, logLevel=acts.logging.VERBOSE
    )

    # Read input space points from input csv files
    evReader = CsvMuonSpacePointReader(
        inputStem="SpacePoints",
        inputDir=os.path.dirname(inDir),
        outputSpacePoints="MuonSpacePoints",
        level=acts.logging.VERBOSE,
    )

    truthReader = CsvMuonSegmentReader(
        inputStem="MuonTruthSegment",
        inputDir=os.path.dirname(inDir),
        outputSegments="MuonTruthSegments",
        level=acts.logging.VERBOSE,
    )

    # add csv reader
    s.addReader(evReader)
    s.addReader(truthReader)
    ### Add the hough seeder algorithm
    seeder = MuonHoughSeeder(
        inSpacePoints=evReader.config.outputSpacePoints,
        inTruthSegments=truthReader.config.outputSegments,
        outHoughMax="MuonHoughSeeds",
        level=acts.logging.VERBOSE,
    )

    s.addAlgorithm(seeder)
    s.run()


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Example script to run ITk seed finding based on CSV space points",
    )
    p.add_argument(
        "--input",
        help="Path to the script's input. By default it's assumed that a ROOT n-tuple is parsed. Otherwise, it's also possible to parse a CSV directory",
    )
    p.add_argument(
        "--isCSV",
        default=False,
        action="store_true",
        help="Flag toggling that the input is a CSV directory",
    )
    p.add_argument("--nEvents", default=100, help="Number of events to run", type=int)

    args = p.parse_args()
    if args.isCSV:
        runHoughFromCsv(inFile=args.input, nEvents=args.nEvents)
    else:
        runHoughFromRoot(args.input, nEvents=args.nEvents)
