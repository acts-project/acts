#!/usr/bin/env python3

"""Build the GeoModel muon mockup and its ACTS tracking geometry."""

import acts
from acts import geomodel as gm
from acts.examples import geomodel as gm_ex


def buildMuonMockup(
    gContext: acts.GeometryContext | None = None,
    input: str = "",
    mockupDetector: str = "Muon",
    stationNames: list[str] | None = None,
    logLevel=acts.logging.INFO,
):
    """Build and return the detector, tracking geometry, field, and mappings."""

    if gContext is None:
        gContext = acts.GeometryContext.dangerouslyDefaultConstruct()
    unit_constants = acts.UnitConstants

    builderConfig = gm_ex.GeoModelMuonMockupBuilder.Config()

    if input:
        geomTree = gm.readFromDb(input)
        builderConfig.stationNames = stationNames or ["BIL", "BML", "BOL"]
    elif mockupDetector == "Muon":
        mockupConfig = gm_ex.GeoMuonMockupExperiment.Config()
        mockupConfig.dumpTree = True
        mockupConfig.dbName = "ActsGeoMS.db"
        mockupConfig.nSectors = 12
        mockupConfig.nEtaStations = 8
        mockupConfig.buildEndcaps = True
        mockupConfig.buildBarrel = True
        mockupBuilder = gm_ex.GeoMuonMockupExperiment(
            mockupConfig, "GeoMockUpMS", logLevel
        )
        builderConfig.stationNames = stationNames or [
            "BI",
            "BM",
            "BO",
            "EAI",
            "EAM",
            "EAO",
            "ECI",
            "ECM",
            "ECO",
        ]
        geomTree = mockupBuilder.constructMS()
    else:
        raise RuntimeError(f"{mockupDetector} not implemented yet")

    factoryConfig = gm.GeoModelDetectorObjectFactory.Config()
    factoryConfig.nameList = [
        "RpcGasGap",
        "MDTDriftGas",
        "TgcGasGap",
        "SmallWheelGasGap",
    ]
    factoryConfig.convertSubVolumes = True
    factoryConfig.convertBox = ["MDT", "RPC", "SmallWheel", "TGC"]

    factoryOptions = gm.GeoModelDetectorObjectFactory.Options()
    factoryOptions.queries = ["Muon"]
    factoryCache = gm.GeoModelDetectorObjectFactory.Cache()
    factory = gm.GeoModelDetectorObjectFactory(factoryConfig, logLevel)
    factory.construct(factoryCache, gContext, geomTree, factoryOptions)
    builderConfig.volumeBoxFPVs = factoryCache.boundingBoxes

    gmDetectorCfg = gm_ex.GeoModelDetector.Config()
    gmDetectorCfg.geoModelTree = geomTree
    gmDetector = gm_ex.GeoModelDetector(gmDetectorCfg)

    trackingGeometryBuilder = gm_ex.GeoModelMuonMockupBuilder(
        builderConfig, "GeoModelMuonMockupBuilder", logLevel
    )

    trackingGeometry = gmDetector.buildTrackingGeometry(gContext, trackingGeometryBuilder)
    return gmDetector, trackingGeometry, factoryConfig.nameList, factoryCache
