import os
import acts
import argparse
from acts import (
    logging,
    GeometryContext,
    CylindricalContainerBuilder,
    DetectorBuilder,
    GeometryIdGenerator,
)

from acts.examples import (
    AlgorithmContext,
    WhiteBoard,
    ObjTrackingGeometryWriter,
)
from acts import geomodel as gm
from acts import examples

from pathlib import Path
from propagation import runPropagation


def runGeant4(
    detector,
    trackingGeometry,
    field,
    outputDir,
    materialMappings=[],
    volumeMappings=[],
    s: acts.examples.Sequencer = None,
    events=100,
    nMuonPerEvt=1,
):
    from acts.examples.simulation import (
        addParticleGun,
        addGeant4,
        EtaConfig,
        MomentumConfig,
        PhiConfig,
        ParticleConfig,
    )
    from pathlib import Path

    s = s or acts.examples.Sequencer(events=events, numThreads=1)
    s.config.logLevel = acts.logging.INFO
    rnd = acts.examples.RandomNumbers()
    u = acts.UnitConstants
    outputDir = Path(outputDir)
    addParticleGun(
        s,
        outputDirRoot=outputDir / "PG",
        momentumConfig=MomentumConfig(50.0 * u.GeV, 60.0 * u.GeV, transverse=True),
        etaConfig=EtaConfig(-1.0, 1.0, uniform=True),
        phiConfig=PhiConfig(0.0, 360.0 * u.degree),
        particleConfig=ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns),
        ),
        multiplicity=nMuonPerEvt,
        rnd=rnd,
        printParticles=False,
    )

    addGeant4(
        s,
        detector,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        outputDirObj=outputDir / "obj",
        rnd=rnd,
        materialMappings=materialMappings,
        volumeMappings=volumeMappings,
    )
    return s


def main():
    from argparse import ArgumentParser

    u = acts.UnitConstants

    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="",  # "/eos/user/c/cimuonsw/GeometryFiles/MockUp.db",
        help="Input SQL file",
    )
    parser.add_argument(
        "--mockupDetector",
        type=str,
        choices=["Muon"],
        help="Predefined mockup detector which is built transiently",
        default="Muon",
    )
    parser.add_argument("--outDir", default="./", help="Output")
    parser.add_argument("--nEvents", default=100, type=int, help="Number of events")

    args = parser.parse_args()

    gContext = acts.GeometryContext()
    logLevel = logging.INFO

    # Create the tracking geometry builder for the muon system
    gmBuilderConfig = gm.GeoModelMuonMockupBuilder.Config()

    # Read the geometry model from the database
    gmTree = None
    ### Use an external geo model file
    if len(args.input):
        gmTree = gm.readFromDb(args.input)
        gmBuilderConfig.stationNames = ["BIL", "BML", "BOL"]

    elif args.mockupDetector == "Muon":
        mockUpCfg = gm.GeoMuonMockupExperiment.Config()
        mockUpCfg.dumpTree = True
        mockUpCfg.dbName = "ActsGeoMS.db"
        mockUpCfg.nSectors = 12
        mockUpCfg.nEtaStations = 8
        mockUpCfg.buildEndcaps = False
        mockUpBuilder = gm.GeoMuonMockupExperiment(mockUpCfg, "GeoMockUpMS", logLevel)
        gmBuilderConfig.stationNames = ["Inner", "Middle", "Outer"]

        gmTree = mockUpBuilder.constructMS()
    else:
        raise RuntimeError(f"{args.mockupDetector} not implemented yet")

    gmFactoryConfig = gm.GeoModelDetectorObjectFactory.Config()
    gmFactoryConfig.nameList = [
        "RpcGasGap",
        "MDTDriftGas",
    ]
    gmFactoryConfig.convertSubVolumes = True
    gmFactoryConfig.convertBox = ["MDT"]

    gmFactory = gm.GeoModelDetectorObjectFactory(gmFactoryConfig, logLevel)
    # The options
    gmFactoryOptions = gm.GeoModelDetectorObjectFactory.Options()
    gmFactoryOptions.queries = ["Muon"]

    # The Cache & construct call
    gmFactoryCache = gm.GeoModelDetectorObjectFactory.Cache()
    gmFactory.construct(gmFactoryCache, gContext, gmTree, gmFactoryOptions)

    gmBuilderConfig.volumeBoxFPVs = gmFactoryCache.boundingBoxes

    gmDetectorCfg = gm.GeoModelDetector.Config()
    gmDetectorCfg.geoModelTree = gmTree
    detector = gm.GeoModelDetector(gmDetectorCfg)

    field = acts.ConstantBField(acts.Vector3(0, 0, 0 * u.T))

    trackingGeometryBuilder = gm.GeoModelMuonMockupBuilder(
        gmBuilderConfig, "GeoModelMuonMockupBuilder", acts.logging.INFO
    )

    trackingGeometry = detector.buildTrackingGeometry(gContext, trackingGeometryBuilder)

    runGeant4(
        detector=detector,
        trackingGeometry=trackingGeometry,
        field=field,
        outputDir=args.outDir,
        volumeMappings=gmFactoryConfig.nameList,
        events=args.nEvents,
    ).run()

    # runPropagation(trackingGeometry, field, args.outDir).run()

    wb = WhiteBoard(acts.logging.INFO)

    context = AlgorithmContext(0, 0, wb, 10)
    obj_dir = Path(args.outDir) / "obj"
    obj_dir.mkdir(exist_ok=True)

    writer = ObjTrackingGeometryWriter(level=acts.logging.INFO, outputDir=str(obj_dir))

    writer.write(context, trackingGeometry)


if __name__ == "__main__":
    main()
