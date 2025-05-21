import acts
import argparse
from acts import (
    logging,
    GeometryContext,
    CylindricalContainerBuilder,
    DetectorBuilder,
    GeometryIdGenerator,
)
from acts import geomodel as gm
from acts import examples


def runGeant4(
    detector,
    trackingGeometry,
    field,
    outputDir,
    materialMappings=[],
    volumeMappings=["MDTDriftGas"],
    s: acts.examples.Sequencer = None,
):
    from acts.examples.simulation import addParticleGun, addGeant4, EtaConfig
    from pathlib import Path

    s = s or acts.examples.Sequencer(events=100, numThreads=1)
    s.config.logLevel = acts.logging.INFO
    rnd = acts.examples.RandomNumbers()
    addParticleGun(
        s,
        EtaConfig(-2.0, 2.0),
        rnd=rnd,
    )
    outputDir = Path(outputDir)
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
        default="/eos/user/c/cimuonsw/GeometryFiles/MockUp.db",
        help="Input SQL file",
    )
    parser.add_argument("--outDir", default="./", help="Output")

    args = parser.parse_args()

    gContext = acts.GeometryContext()
    logLevel = logging.INFO

    # Read the geometry model from the database
    gmTree = gm.readFromDb(args.input)
    gmFactoryConfig = gm.GeoModelDetectorObjectFactory.Config()
    gmFactoryConfig.nameList = [
        "RpcGasGap",
        "MDTDriftGas",
    ]
    gmFactoryConfig.convertSubVolumes = True
    gmFactory = gm.GeoModelDetectorObjectFactory(gmFactoryConfig, logLevel)
    # The options
    gmFactoryOptions = gm.GeoModelDetectorObjectFactory.Options()
    gmFactoryOptions.queries = ["Muon"]
    # The Cache & construct call
    gmFactoryCache = gm.GeoModelDetectorObjectFactory.Cache()
    gmFactory.construct(gmFactoryCache, gContext, gmTree, gmFactoryOptions)

    gmDetectorCfg = gm.GeoModelDetector.Config()
    gmDetectorCfg.geoModelTree = gmTree
    detector = gm.GeoModelDetector(gmDetectorCfg)

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    trackingGeometry = detector.trackingGeometry()
    runGeant4(detector, trackingGeometry, field, args.outDir).run()


if __name__ == "__main__":
    main()
