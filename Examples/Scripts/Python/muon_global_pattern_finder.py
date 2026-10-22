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
):
    # The space points are expressed in their sector frame; the measurement
    # surfaces of the tracking geometry transform them into the global frame
    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    trackingGeometry = TrackingGeometryJsonConverter().fromFile(gctx, geoFile)

    # The per-algorithm timing is only meaningful on a single thread. Building
    # the KDTree of a busy event can underflow in std::nextafter; such floating
    # point exceptions are reported but do not fail the run
    s = acts.examples.Sequencer(
        events=nEvents if nEvents > 0 else None,
        skip=skip,
        numThreads=numThreads,
        failOnUnmaskedFpe=False,
        logLevel=logLevel,
    )

    spReader = RootMuonSpacePointReader(
        filePath=inFile,
        treeName="MuonSpacePoints",
        outputSpacePoints="MuonSpacePoints",
        level=logLevel,
    )
    s.addReader(spReader)

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
    )
