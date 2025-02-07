import sys
from pathlib import Path

import acts.examples
import acts

from acts.examples import (
    CsvSpacePointReader,
    TrackParamsEstimationAlgorithm,
    SeedingPerformanceWriter,
)
from acts.examples.reconstruction import (
    addStandardSeeding,
)

from acts.examples.itk import itkSeedingAlgConfig, InputSpacePointsType


s = acts.examples.Sequencer(
    events=1,
    numThreads=1,
    outputDir = "output"
)

# loggingLevel = acts.logging.INFO
loggingLevel = acts.logging.DEBUG

s.addReader(
    acts.examples.RootAthenaDumpReader(
        level=loggingLevel,
        treename  = "GNN4ITk",
        inputfile = "Dump_GNN4Itk.root",
        onlySpacepoints = True,
        outputPixelSpacePoints = "pixel_spacepoints",
        outputStripSpacePoints = "strip_spacepoints",
        outputSpacePoints = "spacepoints"
    )
)


# run pixel seeding
seeding_pixel = addStandardSeeding(s,
                                   "pixel_spacepoints",
                                   *acts.examples.itk.itkSeedingAlgConfig(
                                       InputSpacePointsType.PixelSpacePoints,
                                       highOccupancyConfig = True),
                                   outputSeeds = "pixel_seeds",
                                   logLevel=loggingLevel
                                   )


# run strip seeding
seeding_strip = addStandardSeeding(s,
                                   "strip_spacepoints",
                                   *acts.examples.itk.itkSeedingAlgConfig(
                                       InputSpacePointsType.StripSpacePoints),
                                   outputSeeds = "strip_seeds",
                                   logLevel=acts.logging.DEBUG
                                   )


s.run()
