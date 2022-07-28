#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse

from acts.examples import (
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
)
from acts.examples.reconstruction import (
    SeedfinderConfigArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
)

import acts

from acts import MaterialMapJsonConverter, UnitConstants as u


def runITk(
    trackingGeometry,
    decorators,
    outputDir: Path,
    events=1,
    outputObj=True,
    outputCsv=False,
    outputJson=False,
):

    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0

        context = AlgorithmContext(ialg, ievt, eventStore)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            csv_dir = outputDir / "csv"
            csv_dir.mkdir(exist_ok=True)
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(csv_dir),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            obj_dir = outputDir / "obj"
            obj_dir.mkdir(exist_ok=True)
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO,
                outputDir=str(obj_dir),
            )
            writer.write(context, trackingGeometry)

        if outputJson:
            json_dir = outputDir / "json"
            json_dir.mkdir(exist_ok=True)
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(json_dir),
                writePerEvent=True,
                writeSensitive=True,
            )
            writer.write(context)

            jmConverterCfg = MaterialMapJsonConverter.Config(
                processSensitives=True,
                processApproaches=True,
                processRepresenting=True,
                processBoundaries=True,
                processVolumes=True,
                processNonMaterial=True,
                context=context.geoContext,
            )

            jmw = JsonMaterialWriter(
                level=acts.logging.VERBOSE,
                converterCfg=jmConverterCfg,
                fileName=str(json_dir / "material-map"),
                writeFormat=JsonFormat.Json,
            )

            jmw.write(trackingGeometry)


def itkSeedingAlgConfig(inputSpacePointsType):

    # variables that do not change for pixel and strip SPs:
    zMax = 3000 * u.mm
    zMin = -3000 * u.mm
    beamPos = (0 * u.mm, 0 * u.mm)
    collisionRegionMin = -200 * u.mm
    collisionRegionMax = 200 * u.mm
    maxSeedsPerSpM = 4
    cotThetaMax = 27.2899
    sigmaScattering = 2
    radLengthPerSeed = 0.1
    minPt = 900 * u.MeV
    bFieldInZ = 2 * u.T
    deltaRMin = 20 * u.mm
    maxPtScattering = float("inf") * u.GeV
    zBinEdges = [
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
    ]  # zBinEdges enables non-equidistant binning in z, in case the binning is not defined the edges are evaluated automatically using equidistant binning
    rRangeMiddleSP = [
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
    ]  # if useVariableMiddleSPRange is set to false, the vector rRangeMiddleSP can be used to define a fixed r range for each z bin: {{rMin, rMax}, ...}. If useVariableMiddleSPRange is set to false and the vector is empty, the cuts won't be applied
    useVariableMiddleSPRange = True  # if useVariableMiddleSPRange is true, the values in rRangeMiddleSP will be calculated based on r values of the SPs and deltaRMiddleSPRange
    binSizeR = 1 * u.mm
    forceRadialSorting = True
    seedConfirmation = True
    centralSeedConfirmationRange = acts.SeedConfirmationRangeConfig(
        zMinSeedConf=-250 * u.mm,
        zMaxSeedConf=250 * u.mm,
        rMaxSeedConf=140 * u.mm,
        nTopForLargeR=1,
        nTopForSmallR=2,
    )  # contains parameters for seed confirmation
    forwardSeedConfirmationRange = acts.SeedConfirmationRangeConfig(
        zMinSeedConf=-3000 * u.mm,
        zMaxSeedConf=3000 * u.mm,
        rMaxSeedConf=140 * u.mm,
        nTopForLargeR=1,
        nTopForSmallR=2,
    )
    compatSeedWeight = 100
    curvatureSortingInFilter = True
    phiMin = 0
    phiMax = 2 * math.pi
    zBinEdges = [
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
    ]  # zBinEdges enables non-equidistant binning in z, in case the binning is not defined the edges are evaluated automatically using equidistant binning
    phiBinDeflectionCoverage = 3
    numPhiNeighbors = 1

    # variables that change for pixel and strip SPs:
    if inputSpacePointsType == "PixelSpacePoints":
        outputSeeds = "PixelSeeds"
        allowSeparateRMax = False
        rMaxGridConfig = 320 * u.mm
        rMaxSeedFinderConfig = rMaxGridConfig
        deltaRMinSP = 6 * u.mm
        deltaRMax = 280 * u.mm
        deltaRMaxTopSP = 280 * u.mm
        deltaRMaxBottomSP = 120 * u.mm
        interactionPointCut = True
        arithmeticAverageCotTheta = False
        deltaZMax = 600 * u.mm
        impactMax = 2 * u.mm
        zBinsCustomLooping = [
            1,
            2,
            3,
            4,
            11,
            10,
            9,
            8,
            6,
            5,
            7,
        ]  # enable custom z looping when searching for SPs, must contain numbers from 1 to the total number of bin in zBinEdges
        skipPreviousTopSP = True
        zBinNeighborsTop = [
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
        ]  # allows to specify the number of neighbors desired for each bin, [-1,1] means one neighbor on the left and one on the right, if the vector is empty the algorithm returns the 8 surrounding bins
        zBinNeighborsBottom = [
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
        ]
        deltaRMiddleMinSPRange = 10 * u.mm
        deltaRMiddleMaxSPRange = 10 * u.mm
        seedConfirmationFilter = True
        impactWeightFactor = 100
        compatSeedLimit = 3
        numSeedIncrement = 10**100  # inf
        seedWeightIncrement = 0
        useDetailedDoubleMeasurementInfo = False
        maxSeedsPerSpMConf = 5
        maxQualitySeedsPerSpMConf = 5
        useDeltaRorTopRadius = True
    else:
        outputSeeds = "StripSeeds"
        allowSeparateRMax = True
        rMaxGridConfig = 1000.0 * u.mm
        rMaxSeedFinderConfig = 1200.0 * u.mm
        deltaRMinSP = 20 * u.mm
        deltaRMax = 600 * u.mm
        deltaRMaxTopSP = 300 * u.mm
        deltaRMaxBottomSP = deltaRMaxTopSP
        interactionPointCut = False
        arithmeticAverageCotTheta = True
        deltaZMax = 900 * u.mm
        impactMax = 20 * u.mm
        zBinsCustomLooping = [6, 7, 5, 8, 4, 9, 3, 10, 2, 11, 1]
        skipPreviousTopSP = False
        zBinNeighborsTop = [
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
        ]
        zBinNeighborsBottom = [
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 2],
            [0, 1],
            [0, 0],
            [-1, 0],
            [-2, 0],
            [-1, 0],
            [-1, 0],
            [-1, 0],
        ]
        deltaRMiddleMinSPRange = 30 * u.mm
        deltaRMiddleMaxSPRange = 150 * u.mm
        seedConfirmationFilter = False
        impactWeightFactor = 1
        compatSeedLimit = 4
        numSeedIncrement = 1
        seedWeightIncrement = 10100
        useDetailedDoubleMeasurementInfo = True
        maxSeedsPerSpMConf = 1000000000
        maxQualitySeedsPerSpMConf = 1000000000
        useDeltaRorTopRadius = False

    # fill namedtuples
    seedfinderConfigArg = SeedfinderConfigArg(
        maxSeedsPerSpM=maxSeedsPerSpM,
        cotThetaMax=cotThetaMax,
        sigmaScattering=sigmaScattering,
        radLengthPerSeed=radLengthPerSeed,
        minPt=minPt,
        bFieldInZ=bFieldInZ,
        impactMax=impactMax,
        interactionPointCut=interactionPointCut,
        arithmeticAverageCotTheta=arithmeticAverageCotTheta,
        deltaZMax=deltaZMax,
        maxPtScattering=maxPtScattering,
        zBinEdges=zBinEdges,
        skipPreviousTopSP=skipPreviousTopSP,
        zBinsCustomLooping=zBinsCustomLooping,
        rRangeMiddleSP=rRangeMiddleSP,
        useVariableMiddleSPRange=useVariableMiddleSPRange,
        binSizeR=binSizeR,
        forceRadialSorting=forceRadialSorting,
        seedConfirmation=seedConfirmation,
        centralSeedConfirmationRange=centralSeedConfirmationRange,
        forwardSeedConfirmationRange=forwardSeedConfirmationRange,
        deltaR=(deltaRMin, deltaRMax),
        deltaRBottomSP=(deltaRMinSP, deltaRMaxBottomSP),
        deltaRTopSP=(deltaRMinSP, deltaRMaxTopSP),
        deltaRMiddleSPRange=(deltaRMiddleMinSPRange, deltaRMiddleMaxSPRange),
        collisionRegion=(collisionRegionMin, collisionRegionMax),
        r=(None, rMaxSeedFinderConfig),
        z=(zMin, zMax),
        beamPos=beamPos,
    )
    seedFilterConfigArg = SeedFilterConfigArg(
        impactWeightFactor=impactWeightFactor,
        compatSeedWeight=compatSeedWeight,
        compatSeedLimit=compatSeedLimit,
        numSeedIncrement=numSeedIncrement,
        seedWeightIncrement=seedWeightIncrement,
        seedConfirmation=seedConfirmation,
        curvatureSortingInFilter=curvatureSortingInFilter,
        maxSeedsPerSpMConf=maxSeedsPerSpMConf,
        maxQualitySeedsPerSpMConf=maxQualitySeedsPerSpMConf,
        useDeltaRorTopRadius=useDeltaRorTopRadius,
    )
    spacePointGridConfigArg = SpacePointGridConfigArg(
        rMax=rMaxGridConfig,
        zBinEdges=zBinEdges,
        phiBinDeflectionCoverage=phiBinDeflectionCoverage,
        phi=(phiMin, phiMax),
        impactMax=impactMax,
    )
    seedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(
        allowSeparateRMax=allowSeparateRMax,
        zBinNeighborsTop=zBinNeighborsTop,
        zBinNeighborsBottom=zBinNeighborsBottom,
        numPhiNeighbors=numPhiNeighbors,
    )

    return (
        seedfinderConfigArg,
        seedFilterConfigArg,
        spacePointGridConfigArg,
        seedingAlgorithmConfigArg,
    )


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Example script to construct the ITk geometry and write it out to CSV and OBJ formats"
    )
    p.add_argument(
        "geo_dir",
        help="Input directory containing the ITk standalone geometry. Get in touch if you don't have this.",
    )
    p.add_argument(
        "--output-dir",
        default=Path.cwd(),
        type=Path,
        help="Directory to write outputs to",
    )
    p.add_argument(
        "--output-csv", action="store_true", help="Write geometry in CSV format."
    )
    p.add_argument(
        "--output-obj", action="store_true", help="Write geometry in OBJ format."
    )
    p.add_argument(
        "--output-json",
        action="store_true",
        help="Write geometry and material in JSON format.",
    )
    p.add_argument(
        "--no-material", action="store_true", help="Decorate material to the geometry"
    )

    args = p.parse_args()
    args.output_dir.mkdir(exist_ok=True, parents=True)

    geo_example_dir = Path(args.geo_dir)
    assert geo_example_dir.exists(), "Detector example input directory missing"
    from acts.examples.itk import buildITkGeometry

    detector, trackingGeometry, decorators = buildITkGeometry(
        geo_example_dir,
        material=not args.no_material,
    )

    runITk(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        outputDir=args.output_dir,
        outputCsv=args.output_csv,
        outputObj=args.output_obj,
        outputJson=args.output_json,
    )
