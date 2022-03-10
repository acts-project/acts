#!/usr/bin/env python3
import os
import argparse

import acts
import acts.examples

from acts.examples import (
     CsvSpacePointReader
)

u = acts.UnitConstants


def runITkSeeding(field, outputDir, s=None):

    # Read input space points from input csv files
    evReader = CsvSpacePointReader(
        level=acts.logging.INFO,
        inputStem="spacepoints",
        inputCollection="pixel",
        inputDir="/acts/Examples/Scripts/Python/CsvSpacePointsOutput_singleMu_100GeV_300evnts/",
        outputSpacePoints="PixelSpacePoints",
    )

    gridConfig = acts.SpacePointGridConfig(
        bFieldInZ=1.997244311 * u.T,
        minPt=900 * u.MeV,
        rMax=320 * u.mm, # pixel: 320 mm, strip: 1000 mm
        zMax=3000 * u.mm,
        zMin=-3000 * u.mm,
        deltaRMax=280 * u.mm, # pixel: 280 mm, strip: 600 mm
        cotThetaMax=27.2899,  # pixel: 27.2899 (4 eta), strip: 900
        impactMax=2 * u.mm, # pixel: 2 mm, strip: 20 mm
        zBinEdges=[-3000., -2500., -1400., -925., -450., -250., 250., 450., 925., 1400., 2500., 3000.], # zBinEdges enables non-equidistant binning in z, in case the binning is not defined the edges are evaluated automatically using equidistant binning
        numPhiNeighbors=1, # number of phiBin neighbors (plus the current bin) that covers the full deflection of a minimum pT particle
    )
    
    seedFinderConfig = acts.SeedfinderConfig(
        rMax=gridConfig.rMax,
        deltaRMin=20 * u.mm,
        deltaRMax=gridConfig.deltaRMax,
        deltaRMinTopSP=6 * u.mm, # pixel: 6 mm, strip: 20 mm
        deltaRMinBottomSP=6 * u.mm, # pixel: 6 mm, strip: 20 mm
        deltaRMaxTopSP=280 * u.mm, # pixel: 280 mm, strip: 3000 mm
        deltaRMaxBottomSP=120 * u.mm, # pixel: 120 mm, strip: 3000 mm
        collisionRegionMin=-200 * u.mm,
        collisionRegionMax=200 * u.mm,
        zMin=gridConfig.zMin,
        zMax=gridConfig.zMax,
        maxSeedsPerSpM=4,
        cotThetaMax=gridConfig.cotThetaMax,
        sigmaScattering=2,
        radLengthPerSeed=0.09804522341059585,
        minPt=gridConfig.minPt,
        bFieldInZ=gridConfig.bFieldInZ,
        beamPos=acts.Vector2(0 * u.mm, 0 * u.mm),
        impactMax=gridConfig.impactMax,
        maxPtScattering=1000000 * u.GeV,
        zBinEdges=gridConfig.zBinEdges,
        enableCutsForSortedSP=True, # enable cotTheta sorting in SeedFinder
        rRangeMiddleSP=[[40., 90.],[40., 200.],[46., 200.],[46., 200.],[46., 250.], [46., 250.], [46., 250.], [46., 200.],[46., 200.], [40., 200.], [40., 90.]], # if useVariableMiddleSPRange is set to false, the vector rRangeMiddleSP can be used to define a fixed r range for each z bin: {{rMin, rMax}, ...}. If useVariableMiddleSPRange is set to false and the vector is empty, the cuts won't be applied
        useVariableMiddleSPRange=True, # if useVariableMiddleSPRange is true, the values in rRangeMiddleSP will be calculated based on r values of the SPs and deltaRMiddleSPRange
        deltaRMiddleSPRange=10,
        seedConfirmation=True,
        centralSeedConfirmationRange=acts.SeedConfirmationRange(250., -250., 140., 1, 2), # contains parameters for seed confirmation (zMinSeedConf, zMaxSeedConf, rMaxSeedConf, nTopForLargeR, nTopForSmallR)
        forwardSeedConfirmationRange=acts.SeedConfirmationRange(3000., -3000., 140., 1, 2),
    )
    
    seedFilterConfig = acts.SeedFilterConfig(
        maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
        deltaRMin=seedFinderConfig.deltaRMin,
        impactWeightFactor = 100,
        compatSeedWeight = 100,
        compatSeedLimit = 3,
    )
    
    seedingAlg = acts.examples.SeedingAlgorithm(
        level=acts.logging.VERBOSE,
        inputSpacePoints=[evReader.config.outputSpacePoints],
        outputSeeds="PixelSeeds",
        outputProtoTracks="prototracks",
        gridConfig=gridConfig,
        seedFinderConfig=seedFinderConfig,
        seedFilterConfig=seedFilterConfig,
        zBinNeighborsTop=[[0, 0],[-1, 0],[-1, 0],[-1, 0],[-1, 0],[-1, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 0]], # allows to specify the number of neighbors desired for each bin, [-1,1] means one neighbor on the left and one on the right, if the vector is empty the algorithm returns the 8 surrounding bins
        zBinNeighborsBottom=[[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 0],[-1, 0],[-1, 0],[-1, 0],[-1, 0],[-1, 0]],
    )

    s = s or acts.examples.Sequencer(
        events=1, numThreads=-1, logLevel=acts.logging.INFO
    )

    s.addReader(evReader)
    s.addAlgorithm(seedingAlg)

    return s


if "__main__" == __name__:

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runITkSeeding(field, outputDir=os.getcwd()).run()
