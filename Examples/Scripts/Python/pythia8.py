#!/usr/bin/env python3
from pathlib import Path

import acts
import acts.examples
from acts.examples.simulation import addPythia8

u = acts.UnitConstants


def runPythia8(
    outputDir,
    outputRoot: bool = True,
    outputCsv: bool = True,
    s: acts.examples.Sequencer = None,
    vtxGen=None,
):
    # Preliminaries
    rnd = acts.examples.RandomNumbers()
    outputDir = Path(outputDir)

    # Sequencer
    s = s or acts.examples.Sequencer(
        events=10, numThreads=-1, logLevel=acts.logging.INFO
    )

    addPythia8(
        s,
        rnd=rnd,
        outputDirCsv=outputDir / "csv" if outputCsv else None,
        outputDirRoot=outputDir if outputRoot else None,
        vtxGen=vtxGen,
    )

    return s


if "__main__" == __name__:
    runPythia8(Path.cwd()).run()
