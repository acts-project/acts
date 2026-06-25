import acts
import argparse

from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants


def main():

    import acts

    global __geant4Handle
    global __pmr

    p = argparse.ArgumentParser()

    p.add_argument(
        "-i", "--input", type=str, default="", help="Input file(s) for the geometry"
    )

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
        type=lambda x: acts.logging.Level.__members__[x.upper()],
        default=acts.logging.INFO,
        metavar="|".join(acts.logging.Level.__members__.keys()),
        help="Log level for the geometry construction and conversion",
    )

    # Output the detector to json file, default is false
    p.add_argument(
        "--output-json",
        action="store_true",
        help="Output the detector to a json file for visualization",
    )

    # Run the consistency check for the detray geometry, default is false
    p.add_argument(
        "--detray-consistency-check",
        action="store_true",
        help="Run the consistency check for the detray geometry",
    )

    args = p.parse_args()

    prfx = args.output + "_" if args.output != "" else ""

    gContext = acts.GeometryContext.dangerouslyDefaultConstruct()

    logLevel = args.log_level

    # Material decoration for reconstruction geometry
    materialDecorator = None
    if args.map != "":
        print(">>> Loading a material decorator from file:", args.map)
        materialDecorator = acts.IMaterialDecorator.fromFile(args.map)

    trackingGeometry = None
    detectorStore = {}

    # Build from ODD unless an input file is provided for a mode that supports it
    build_from_odd = args.input == "" or args.geo_mode in ("gen1", "geant4")

    if build_from_odd:
        buildGen3 = args.geo_mode == "gen3" or args.geo_mode == "detray"
        with getOpenDataDetector(gen3=buildGen3, logLevel=logLevel) as detector:
            trackingGeometry = detector.trackingGeometry()
            detectorStore["Detector"] = detector
            detectorStore["Volume"] = trackingGeometry.highestTrackingVolume
            detectorStore["SurfaceByIdentifier"] = trackingGeometry.geoIdSurfaceMap()

    print(">>> Test mode is :", args.geo_mode)

    if args.geo_mode == "gen3":
        from pathlib import Path
        from acts.json import TrackingGeometryJsonConverter

        converter = TrackingGeometryJsonConverter(level=logLevel)

        if args.input != "":
            print(">>> Reading tracking geometry from", args.input)
            trackingGeometry = converter.fromJson(
                gContext, Path(args.input).read_text()
            )
            print(">>> Read tracking geometry from", args.input)

        if args.output_json:
            print(">>> Outputting the tracking geometry to json file ...")
            json_str = converter.toJson(gContext, trackingGeometry)
            out_path = prfx + "tracking-geometry.json"
            Path(out_path).write_text(json_str)
            print(">>> Written to", out_path)

    # check if the mode does not contain geant4
    elif args.geo_mode == "detray":
        import glob
        import acts.vecmem, acts.detray
        import acts.examples.detray

        __pmr = acts.vecmem.HostMemoryResource()

        if args.input != "":
            files = glob.glob(args.input.rstrip("/") + "/*.json")
            print(">>> Reading detray geometry from", args.input, "->", files)
            detrayDetector, detrayNames = acts.detray.readODD(__pmr, files)
        else:
            payloadConfig = acts.detray.DetrayPayloadConverter.Config()
            payloadConfig.beampipeVolume = trackingGeometry.findVolumeByName(
                "BeamPipe"
            )
            payloadConverter = acts.detray.DetrayPayloadConverter(
                payloadConfig, logLevel
            )

            converterConfig = acts.detray.DetrayGeometryConverter.Config()
            converterConfig.payloadConverter = payloadConverter
            converter = acts.detray.DetrayGeometryConverter(
                converterConfig, logLevel
            )

            detrayGeometry = converter.convert(
                __pmr, gContext, trackingGeometry, detectorName="odd"
            )
            detrayDetector, detrayNames = (
                detrayGeometry.detector,
                detrayGeometry.names,
            )

        if args.detray_consistency_check:
            detrayDetector.checkConsistency()

        if args.output_json:
            print(">>> Outputting the detray geometry to json file ...")
            detray_out = prfx + "detray/"
            detrayDetector.writeToJson(detrayNames, detray_out)
            print(">>> Written to", detray_out)

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
