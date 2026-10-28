#!/usr/bin/env python3
"""Simulate event ranges in independent Geant4 processes, with CSV/ROOT output."""

import argparse
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory

from acts.examples.geant4_multiprocessing import (
    mergeGeant4RootFiles,
    runGeant4Processes,
)


def runGeant4EventRange(
    beginEvent,
    endEvent,
    *,
    outputDir,
    rootDir=None,
    detectorFactory=None,
    seed=42,
    etaRange=(-2.0, 2.0),
):
    import acts
    import acts.examples
    from acts.examples.simulation import addParticleGun, addGeant4, EtaConfig

    if detectorFactory is None:
        from acts.examples.odd import getOpenDataDetector

        detectorFactory = getOpenDataDetector
    # Construct all geometry and Geant4 state in this fresh process.
    detector = detectorFactory()
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))
    rnd = acts.examples.RandomNumbers(seed=seed)
    s = acts.examples.Sequencer(
        events=endEvent - beginEvent, skip=beginEvent, numThreads=1
    )
    for decorator in detector.contextDecorators():
        s.addContextDecorator(decorator)
    addParticleGun(
        s,
        EtaConfig(*etaRange),
        rnd=rnd,
        outputDirCsv=Path(outputDir) / "csv",
    )
    addGeant4(
        s,
        detector,
        trackingGeometry,
        field,
        outputDirCsv=Path(outputDir) / "csv",
        outputDirRoot=(
            Path(rootDir) / f"events_{beginEvent}_{endEvent}"
            if rootDir is not None
            else None
        ),
        rnd=rnd,
    )
    s.run()


def runGeant4Parallel(
    *,
    events=100,
    processes=8,
    skip=0,
    outputDir=Path.cwd(),
    outputRoot=False,
    detectorFactory=None,
    seed=42,
    etaRange=(-2.0, 2.0),
):
    """Create one detector/run per process; retain global event IDs and seeds."""
    outputDir = Path(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)
    (outputDir / "csv").mkdir(exist_ok=True)
    with TemporaryDirectory(prefix="geant4-workers-", dir=outputDir) as temporary:
        ranges = runGeant4Processes(
            partial(
                runGeant4EventRange,
                outputDir=outputDir,
                rootDir=temporary if outputRoot else None,
                detectorFactory=detectorFactory,
                seed=seed,
                etaRange=etaRange,
            ),
            events=events,
            processes=processes,
            skip=skip,
        )
        if outputRoot:
            names = ("particles_simulation.root", "hits.root")
            for name in names:
                mergeGeant4RootFiles(
                    [
                        Path(temporary) / f"events_{begin}_{end}" / name
                        for begin, end in ranges
                    ],
                    Path(temporary) / name,
                )
            for name in names:
                (Path(temporary) / name).replace(outputDir / name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-n", "--events", type=int, default=100)
    parser.add_argument("-j", "--processes", type=int, default=8)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("-o", "--output", type=Path, default=Path.cwd())
    parser.add_argument(
        "--root", action="store_true", help="Also write merged ROOT output"
    )
    args = parser.parse_args()
    runGeant4Parallel(
        events=args.events,
        processes=args.processes,
        skip=args.skip,
        seed=args.seed,
        outputDir=args.output,
        outputRoot=args.root,
    )


if __name__ == "__main__":
    main()
