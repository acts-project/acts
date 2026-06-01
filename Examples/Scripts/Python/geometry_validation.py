import acts
import argparse

from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants


def main():

    import acts

    global __geant4Handle
    global __pmr

    p = argparse.ArgumentParser()

    p.add_argument("-i", "--input", type=str, default="", help="Input SQL file")

    p.add_argument("-o", "--output", type=str, default="", help="Output prefix")

    p.add_argument(
        "-m", "--map", type=str, default="", help="Input file for the material map"
    )

    # The geometry modes are
    # gen1: Gen1 detector with Gen1 navigator and propagator
    # gen3: Gen3 detector with Gen3 navigator and propagator
    # detray: Gen3->detray detector with detray navigator and propagator
    # geant4: Geant4 navigator and propagator with gen3 surface matching
    p.add_argument(
        "--geo-mode",
        type=str,
        default="gen3",
        choices=["gen1", "gen3", "detray", "geant4"],
        help="Convert to detray detector and run detray navigation and propagation",
    )

    # Set log level for the geometry construction and conversion
    p.add_argument(
        "--log-level",
        type=int,
        default=acts.logging.INFO,
        help="Log level for the geometry construction and conversion",
    )

    args = p.parse_args()

    prfx = args.output + "_" if args.output != "" else ""

    gContext = acts.GeometryContext.dangerouslyDefaultConstruct()

    # translate integet log level to acts logging level
    logLevel = acts.examples.getLogLevel(args.log_level)

    # Material decoration for reconstruction geometry
    materialDecorator = None
    if args.map != "":
        print(">>> Loading a material decorator from file:", args.map)
        materialDecorator = acts.IMaterialDecorator.fromFile(args.map)

    trackingGeometry = None
    detectorStore = {}

    buildGen3 = args.geo_mode == "gen3" or args.geo_mode == "detray"

    with getOpenDataDetector(gen3=buildGen3, logLevel=logLevel) as detector:
        trackingGeometry = detector.trackingGeometry()
        detectorStore["Detector"] = detector
        detectorStore["Volume"] = trackingGeometry.highestTrackingVolume
        detectorStore["SurfaceByIdentifier"] = trackingGeometry.geoIdSurfaceMap()

    print(">>> Test mode is :", args.geo_mode)
    # check if the mode does not contain geant4
    if args.geo_mode == "detray":
        import acts.vecmem, acts.detray
        import acts.examples.detray

        __pmr = acts.vecmem.HostMemoryResource()

        detrayGeometry = acts.detray.convertODD(
            __pmr,
            gContext,
            trackingGeometry,
            beampipeVolumeName="BeamPipe",
            logLevel=logLevel,
        )
    elif args.geo_mode == "geant4":

        print(">>> Building Geant4 detector ...")
        detector = detectorStore["Detector"]

        import acts.examples.geant4

        from acts.examples.geant4 import (
            Geant4Simulation,
            Geant4ConstructionOptions,
            SensitiveSurfaceMapper,
        )

        smmConfig = SensitiveSurfaceMapper.Config()
        smmConfig.volumeMappings = []
        smmConfig.materialMappings = ["Silicon"]
        sensitiveMapper = SensitiveSurfaceMapper.create(
            smmConfig, acts.logging.INFO, trackingGeometry
        )

        detectorConstructionOptions = acts.examples.geant4.Geant4ConstructionOptions()

        # Specify the physics list, kill volume and other options for the Geant4 simulation
        physicsList = "FTFP_BERT"  # "MaterialPhysicsList"
        killVolume = detectorStore["Volume"]
        killAfterTime = float("inf")
        bfield = None

        inputParticles = "particles_generated"
        outputParticles = "particles_final"


if "__main__" == __name__:
    main()
