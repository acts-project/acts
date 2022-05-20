#!/usr/bin/env python3
import os
import argparse
import tempfile

import acts
import acts.examples

from acts.examples import CsvSpacePointReader

u = acts.UnitConstants


def runITkSeeding(field, csvInputDir, outputDir, s=None):

    # Read input space points from input csv files
    evReader = CsvSpacePointReader(
        level=acts.logging.INFO,
        inputStem="spacepoints",
        inputCollection="pixel",
        inputDir=csvInputDir,
        outputSpacePoints="PixelSpacePoints",
        extendCollection=False,
    )

    gridConfig = acts.SpacePointGridConfig(
        bFieldInZ=2 * u.T,
        minPt=900 * u.MeV,
        rMax=320 * u.mm,  # pixel: 320 mm, strip: 1000 mm
        zMax=3000 * u.mm,
        zMin=-3000 * u.mm,
        deltaRMax=280 * u.mm,  # pixel: 280 mm, strip: 600 mm
        cotThetaMax=27.2899,  # pixel: 27.2899 (4 eta)
        impactMax=2 * u.mm,  # pixel: 2 mm, strip: 20 mm
        zBinEdges=[
            -3000.0,
            -2500.0,
            -1400.0,
            -925.0,
            -450.0,
            -250.0,
            250.0,
            450.0,
            925.0,
            1400.0,
            2500.0,
            3000.0,
        ],  # zBinEdges enables non-equidistant binning in z, in case the binning is not defined the edges are evaluated automatically using equidistant binning
        phiBinDeflectionCoverage=3,
    )

    seedFinderConfig = acts.SeedfinderConfig(
        rMax=gridConfig.rMax,
        deltaRMin=20 * u.mm,
        deltaRMax=gridConfig.deltaRMax,
        deltaRMinTopSP=6 * u.mm,  # pixel: 6 mm, strip: 20 mm
        deltaRMinBottomSP=6 * u.mm,  # pixel: 6 mm, strip: 20 mm
        deltaRMaxTopSP=280 * u.mm,  # pixel: 280 mm, strip: 3000 mm
        deltaRMaxBottomSP=120 * u.mm,  # pixel: 120 mm, strip: 3000 mm
        collisionRegionMin=-200 * u.mm,
        collisionRegionMax=200 * u.mm,
        zMin=gridConfig.zMin,
        zMax=gridConfig.zMax,
        maxSeedsPerSpM=4,
        cotThetaMax=gridConfig.cotThetaMax,
        sigmaScattering=2,
        radLengthPerSeed=0.1,
        minPt=gridConfig.minPt,
        bFieldInZ=gridConfig.bFieldInZ,
        beamPos=acts.Vector2(0 * u.mm, 0 * u.mm),
        impactMax=gridConfig.impactMax,
        maxPtScattering=float("inf") * u.GeV,
        deltaZMax=900 * u.mm,
        interactionPointCut=True,
        zBinEdges=gridConfig.zBinEdges,
        skipPreviousTopSP=True,
        rRangeMiddleSP=[
            [40.0, 90.0],
            [40.0, 200.0],
            [46.0, 200.0],
            [46.0, 200.0],
            [46.0, 250.0],
            [46.0, 250.0],
            [46.0, 250.0],
            [46.0, 200.0],
            [46.0, 200.0],
            [40.0, 200.0],
            [40.0, 90.0],
        ],  # if useVariableMiddleSPRange is set to false, the vector rRangeMiddleSP can be used to define a fixed r range for each z bin: {{rMin, rMax}, ...}. If useVariableMiddleSPRange is set to false and the vector is empty, the cuts won't be applied
        useVariableMiddleSPRange=True,  # if useVariableMiddleSPRange is true, the values in rRangeMiddleSP will be calculated based on r values of the SPs and deltaRMiddleSPRange
        deltaRMiddleMinSPRange=10,
        deltaRMiddleMaxSPRange=10,
        seedConfirmation=True,
        centralSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=250 * u.mm,
            zMaxSeedConf=250 * u.mm,
            rMaxSeedConf=140 * u.mm,
            nTopForLargeR=1,
            nTopForSmallR=2,
        ),  # contains parameters for seed confirmation
        forwardSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=3000 * u.mm,
            zMaxSeedConf=-3000 * u.mm,
            rMaxSeedConf=140 * u.mm,
            nTopForLargeR=1,
            nTopForSmallR=2,
        ),
        useDetailedDoubleMeasurementInfo=False,
    )

    seedFilterConfig = acts.SeedFilterConfig(
        maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
        deltaRMin=seedFinderConfig.deltaRMin,
        impactWeightFactor=100,
        compatSeedWeight=100,
        compatSeedLimit=3,
        curvatureSortingInFilter=True,
        seedConfirmation=True,
        centralSeedConfirmationRange=seedFinderConfig.centralSeedConfirmationRange,
        forwardSeedConfirmationRange=seedFinderConfig.forwardSeedConfirmationRange,
        useDeltaRorTopRadius=True,
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=acts.logging.VERBOSE,
        inputSpacePoints=[evReader.config.outputSpacePoints],
        outputSeeds="PixelSeeds",
        outputProtoTracks="prototracks",
        allowSeparateRMax=False,
        gridConfig=gridConfig,
        seedFinderConfig=seedFinderConfig,
        seedFilterConfig=seedFilterConfig,
        zBinNeighborsTop=[
            [0, 0],
            [-1, 0],
            [-1, 0],
            [-1, 0],
            [-1, 0],
            [-1, 1],
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 0],
        ],  # allows to specify the number of neighbors desired for each bin, [-1,1] means one neighbor on the left and one on the right, if the vector is empty the algorithm returns the 8 surrounding bins
        zBinNeighborsBottom=[
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 0],
            [-1, 0],
            [-1, 0],
            [-1, 0],
            [-1, 0],
            [-1, 0],
        ],
        numPhiNeighbors=1,
    )

    s = s or acts.examples.Sequencer(
        events=1, numThreads=-1, logLevel=acts.logging.INFO
    )

    s.addReader(evReader)
    s.addAlgorithm(seedingAlg)

    return s


if "__main__" == __name__:

    # create temporary file
    with tempfile.TemporaryDirectory() as tmpdirname:
        temp = open(tmpdirname + "/event000000000-spacepoints_pixel.csv", "w+t")
        print(
            "created temporary file: "
            + tmpdirname
            + "/event000000000-spacepoints_pixel.csv"
        )
        temp.write(
            "measurement_id,sp_type,module_idhash,sp_x,sp_y,sp_z,sp_radius,sp_covr,sp_covz,sp_topHalfStripLength,sp_bottomHalfStripLength,sp_topStripDirection[0],sp_topStripDirection[1],sp_topStripDirection[2],sp_bottomStripDirection[0],sp_bottomStripDirection[1],sp_bottomStripDirection[2],sp_stripCenterDistance[0],sp_stripCenterDistance[1],sp_stripCenterDistance[2],sp_bottomStripCenterPosition[0],sp_bottomStripCenterPosition[1],sp_bottomStripCenterPosition[2]\n 1,0,3139,32.67557144165039,-5.311902523040771,-47.65000152587891,33.10452270507812,0.05999999865889549,0.02999880164861679,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n 2,0,3422,95.14442443847656,-15.46361255645752,-52.125,96.39286804199219,0.05999999865889549,0.01687432639300823,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n 3,0,3650,102.8257064819336,-16.71612739562988,-52.67499923706055,104.1755981445312,0.05999999865889549,0.001875000074505806,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n 4,0,4223,159.4266204833984,-25.91166687011719,-56.75,161.5186157226562,0.05999999865889549,0.02999880164861679,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n 5,0,5015,224.07958984375,-36.37123107910156,-61.40000152587891,227.0121765136719,0.05999999865889549,0.007499700412154198,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n 6,0,6023,284.1485595703125,-46.0638542175293,-65.72499847412109,287.8580932617188,0.05999999865889549,0.001875000074505806,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
        )
        temp.read()

        # set magnetic field
        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

        # run seeding
        runITkSeeding(field, os.path.dirname(temp.name), outputDir=os.getcwd()).run()
