#!/usr/bin/env python3
import argparse

import acts
import acts.examples
from acts.examples.root import RootMuonSpacePointReader
from acts.json import TrackingGeometryJsonConverter


def runGlobalPatternFinder(
    inFile: str,
    geoFile: str,
    nEvents: int,
    logLevel: acts.logging.Level,
    numThreads: int = 1,
    skip: int = 0,
    timingDir: str = "",
    timingFile: str = "timing.csv",
):
    # The muon space points are stored in the frame of their spectrometer sector.
    # Reaching the global frame needs the measurement surfaces, which are read
    # from the tracking geometry JSON written alongside the n-tuple.
    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    trackingGeometry = TrackingGeometryJsonConverter().fromFile(gctx, geoFile)

    # The sequencer drives the event loop: for every event it calls each
    # reader, then each algorithm, in the order they were added. Data is
    # passed between them through the event store under string keys.
    #
    # It writes the per algorithm wall time to timingFile on its own. That
    # measurement only means anything on a single thread, which is why one
    # thread is the default: more of them measure the throughput of the machine
    # instead of the cost of the algorithm.
    s = acts.examples.Sequencer(
        events=nEvents if nEvents > 0 else None,
        skip=skip,
        numThreads=numThreads,
        outputDir=timingDir,
        outputTimingFile=timingFile,
        logLevel=logLevel,
    )

    # Reads the space points and writes them into the event store as a
    # MuonSpacePointContainer under the key given by outputSpacePoints. The tree
    # holding them is called "MuonSpacePoints", which differs from the reader's
    # default.
    spReader = RootMuonSpacePointReader(
        filePath=inFile,
        treeName="MuonSpacePoints",
        outputSpacePoints="MuonSpacePoints",
        level=logLevel,
    )
    s.addReader(spReader)

    # Reads the space points back via inSpacePoints and writes the found
    # patterns under outPatterns. Every Config field not set here keeps its
    # default from the C++ header.
    patternFinder = acts.examples.GlobalPatternFinderAlgorithm(
        inSpacePoints=spReader.config.outputSpacePoints,
        outPatterns="MuonGlobalPatterns",
        trackingGeometry=trackingGeometry,
        level=logLevel,
    )
    s.addAlgorithm(patternFinder)

    s.run()


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Run the muon global pattern finder on muon space points",
    )
    p.add_argument(
        "--input",
        required=True,
        help="Path to the ROOT n-tuple with the MuonSpacePoints tree",
    )
    p.add_argument(
        "--geometry",
        required=True,
        help="Path to the tracking geometry JSON matching the n-tuple",
    )
    p.add_argument(
        "--nEvents",
        default=100,
        type=int,
        help="Number of events to run, 0 for all of them",
    )
    p.add_argument(
        "--threads",
        default=1,
        type=int,
        help="Parallel threads, -1 for one per core. Leave at 1 to time the "
        "algorithm rather than the machine",
    )
    p.add_argument(
        "--skip", default=0, type=int, help="Events to skip before the first one"
    )
    p.add_argument(
        "--timingDir",
        default="",
        help="Directory of the timing file, empty for the working directory",
    )
    p.add_argument(
        "--timingFile",
        default="timing.csv",
        help="Name of the file the per algorithm wall time is written to",
    )
    p.add_argument(
        "--logLevel",
        default="INFO",
        choices=["VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR"],
        help="Log level of the sequencer, reader and algorithm",
    )

    args = p.parse_args()
    runGlobalPatternFinder(
        inFile=args.input,
        geoFile=args.geometry,
        nEvents=args.nEvents,
        logLevel=getattr(acts.logging, args.logLevel),
        numThreads=args.threads,
        skip=args.skip,
        timingDir=args.timingDir,
        timingFile=args.timingFile,
    )
