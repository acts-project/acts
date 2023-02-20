#!/usr/bin/env python3
from pathlib import Path
from multiprocessing import Pool
from functools import partial

import numpy as np

# This script runs a Geant4 simulation in parallel by passing chunks of events to subprocesses.
# This is a workaround to achieve parallel processing even though Geant4 is not thread-save
# and thus the internal parallelism of the ACTS examples framework cannot be used.
#
# Note that:
# * This should give equivalent results to a sequential run if the RNG is initialized with
#   the same seed in all runs
# * This works only for csv outputs, since they write the results in one file per event
#   (the naming of the csv-files should be equivalent to a sequential run)


def runGeant4EventRange(events, outputDir):
    import acts
    import acts.examples
    from acts.examples.simulation import addParticleGun, addGeant4, EtaConfig
    from acts.examples.odd import getOpenDataDetector
    from common import getOpenDataDetectorDirectory

    u = acts.UnitConstants

    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory()
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(events=len(events), skip=events[0], numThreads=1)

    outputDir = Path(outputDir)
    addParticleGun(
        s,
        EtaConfig(-2.0, 2.0),
        rnd=rnd,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=None,
    )
    addGeant4(
        s,
        detector,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=None,
        rnd=rnd,
    )

    s.run()
    del s


if "__main__" == __name__:
    n_events = 100
    n_jobs = 8

    outputDir = Path.cwd()

    with Pool(n_jobs) as p:
        p.map(
            partial(runGeant4EventRange, outputDir=outputDir),
            np.array_split(np.arange(n_events), n_jobs),
        )
