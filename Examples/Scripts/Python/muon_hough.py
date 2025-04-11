#!/usr/bin/env python3
import os
import argparse

import acts

from acts.examples import CsvMuonSpacePointReader, CsvMuonSegmentReader, MuonHoughSeeder

# from acts.examples.reconstruction import (
#     addSeeding,
#     addStandardSeeding,
#     SeedingAlgorithm,
# )

# from acts.examples.itk import itkSeedingAlgConfig, InputSpacePointsType

u = acts.UnitConstants
rnd = acts.examples.RandomNumbers(seed=42)


def runHoughFromCsv(inDir):
    # create temporary file with pixel SPs and run the seeding

    s = acts.examples.Sequencer(events=8, numThreads=1, logLevel=acts.logging.VERBOSE)

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
        description="Example script to run ITk seed finding based on CSV spacepoints",
    )
    p.add_argument(
        "indir",
        help="Input directory containing the ITk standalone geometry. Get in touch if you don't have this.",
    )

    args = p.parse_args()

    runHoughFromCsv(args.indir)
