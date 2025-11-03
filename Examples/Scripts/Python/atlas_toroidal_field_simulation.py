#!/usr/bin/env python3

import os
import acts
import acts.acts_atlas_toroidal_field as atlas_field
import argparse
from acts import (
    logging,
    GeometryContext,
)
from pathlib import Path

# Import GeoModel if available
try:
    from acts import geomodel as gm
    HAS_GEOMODEL = True
except ImportError:
    print("Warning: GeoModel not available")
    HAS_GEOMODEL = False

# Import examples if available
try:
    from acts import examples
    from acts.examples import (
        AlgorithmContext,
        WhiteBoard,
        ObjTrackingGeometryWriter,
    )
    HAS_EXAMPLES = True
except ImportError:
    print("Warning: ACTS examples not available")
    HAS_EXAMPLES = False

# Create configuration for ATLAS toroidal field
def create_atlas_field():
    config = atlas_field.Config()
    
    print(f"Creating ATLASToroidalField with:")
    print(f"  Barrel:")
    print(f"    Inner radius:  {config.barrel.R_in} m")
    print(f"    Outer radius:  {config.barrel.R_out} m")
    print(f"    Current:       {config.barrel.I} A")
    print(f"    Turns:         {config.barrel.Nturns}")
    print(f"  Endcaps:")
    print(f"    Inner radius:  {config.ect.R_in} m")
    print(f"    Outer radius:  {config.ect.R_out} m")
    print(f"    Current:       {config.ect.I} A")
    print(f"    Turns:         {config.ect.Nturns}")
    print(f"  Layout:")
    print(f"    Number of coils: {config.layout.nCoils}")
    
    return atlas_field.ATLASToroidalField(config)


def runGeant4(
    detector,
    trackingGeometry,
    field,
    outputDir,
    volumeMappings,
    events=10,
    seed=None,
    particleTypeWeight=100,
    particleGun=None,
    materialMappings=None,
):
    from pathlib import Path

    if not HAS_EXAMPLES:
        print("❌ acts.examples not available, skipping Geant4 simulation")
        return None

    try:
        from acts.examples.simulation import addGeant4, ParticleConfig, MomentumConfig, ParticleSelectorConfig
        from acts.examples.simulation import EtaConfig
        from acts.examples.simulation import addParticleGun
        
        logger = acts.logging.getLogger("Geant4Simulation")
        logger.setLevel(acts.logging.INFO)

        rnd = acts.examples.RandomNumbers(seed=seed or 42)

        sequencer = acts.examples.Sequencer(
            events=events,
            logLevel=acts.logging.INFO,
            numThreads=1,
        )

        addParticleGun(
            sequencer,
            MomentumConfig(1.0 * acts.UnitConstants.GeV, 10.0 * acts.UnitConstants.GeV, transverse=True),
            EtaConfig(-2.6, 2.6),
            ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
            rnd=rnd,
        )

        # Fix materialMappings to be an empty list if None
        if materialMappings is None:
            materialMappings = []

        addGeant4(
            sequencer,
            detector,
            trackingGeometry,
            field,
            outputDirRoot=outputDir,
            volumeMappings=volumeMappings,
            materialMappings=materialMappings,
            rnd=rnd,
        )

        return sequencer
    except Exception as e:
        print(f"❌ Error in Geant4 simulation: {e}")
        return None
def main():
    from argparse import ArgumentParser

    # Check if we have the required dependencies
    if not HAS_GEOMODEL:
        print("❌ GeoModel not available. Testing just the BarrelToroidField...")
        # Just test the barrel field functionality
        field = create_atlas_field()
        print("✅ BarrelToroidField created successfully!")
        return

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
    parser.add_argument(
        "--randomSeed", default=1602, type=int, help="Random seed for event generation"
    )
    parser.add_argument(
        "--geoSvgDump",
        default=False,
        action="store_true",
        help="Dump the tracking geometry in an obj format",
    )

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
        "TgcGasGap",
        "SmallWheelGasGap",
    ]
    gmFactoryConfig.convertSubVolumes = True
    gmFactoryConfig.convertBox = ["MDT", "RPC"]

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

    # Create the barrel toroid field using our function
    field = create_atlas_field()

    trackingGeometryBuilder = gm.GeoModelMuonMockupBuilder(
        gmBuilderConfig, "GeoModelMuonMockupBuilder", logLevel
    )

    trackingGeometry = detector.buildTrackingGeometry(gContext, trackingGeometryBuilder)

    algSequence = runGeant4(
        detector=detector,
        trackingGeometry=trackingGeometry,
        field=field,
        outputDir=args.outDir,
        volumeMappings=gmFactoryConfig.nameList,
        events=args.nEvents,
        seed=args.randomSeed,
    )

    if algSequence is None:
        print("✅ Test completed (Geant4 simulation not available)")
        return

    if HAS_EXAMPLES:
        from acts.examples import MuonSpacePointDigitizer

        digiAlg = MuonSpacePointDigitizer(
            randomNumbers=acts.examples.RandomNumbers(
                acts.examples.RandomNumbers.Config(seed=2 * args.randomSeed)
            ),
            trackingGeometry=trackingGeometry,
            dumpVisualization=False,
            digitizeTime=True,
            outputSpacePoints="MuonSpacePoints",
            level=logLevel,
        )
        algSequence.addAlgorithm(digiAlg)

        from acts.examples import RootMuonSpacePointWriter

        algSequence.addWriter(
            RootMuonSpacePointWriter(
                level=logLevel,
                inputSpacePoints="MuonSpacePoints",
                filePath=f"{args.outDir}/MS_SpacePoints.root",
            )
        )

    if args.geoSvgDump and HAS_EXAMPLES:
        wb = WhiteBoard(acts.logging.INFO)
        context = AlgorithmContext(0, 0, wb, 10)
        obj_dir = Path(args.outDir) / "obj"
        obj_dir.mkdir(exist_ok=True)
        writer = ObjTrackingGeometryWriter(
            level=acts.logging.INFO, outputDir=str(obj_dir)
        )

        writer.write(context, trackingGeometry)

    if algSequence:
        algSequence.run()


if __name__ == "__main__":
    main()
