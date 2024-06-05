import acts
import argparse
from acts import logging, GeometryContext
from acts import geomodel as gm
from acts import examples 

def main():

    p = argparse.ArgumentParser()

    p.add_argument(
        "-i", "--input", type=str, default="", help="Input SQL file"
    )

    p.add_argument(
        "-q", "--queries", type=str, nargs="+", default= "GeoModelXML", help="List of Queries for Published full phys volumes"
    )

    p.add_argument(
    "--output-obj",
    help="Write the surfaces to OBJ files",
    action="store_true",
    default=False,
    )

    args = p.parse_args()

    gContext = acts.GeometryContext()

    # Read the geometry model from the database
    gmTree = acts.geomodel.readFromDb(args.input)

    gmFactoryConfig = gm.GeoModelDetectorSurfaceFactory.Config()
    gmFactoryConfig.shapeConverters = [
        gm.GeoBoxConverter(),
        gm.GeoTrdConverter(),
        gm.GeoIntersectionAnnulusConverter(),
        gm.GeoShiftConverter(),
        gm.GeoUnionDoubleTrdConverter(),
    ]
    gmFactoryConfig.nameList = [
        # "Sensor",
        # "FwdSen",
        # "Wheel",
        # "Side",
        # "Pixel",
        # "Module", # Catches a lot
        # "ModuleBrl",
        # "FwdSensor_Side#1_2_0_0_33",
        # "ECSensor0",
    ]
    gmFactoryConfig.materialList = [
        "std::Silicon",
    ]

    gmFactory = gm.GeoModelDetectorSurfaceFactory(gmFactoryConfig, logging.VERBOSE)
    # The options
    gmFactoryOptions = gm.GeoModelDetectorSurfaceFactory.Options()
    gmFactoryOptions.queries = args.queries
    # The Cache & construct call
    gmFactoryCache = gm.GeoModelDetectorSurfaceFactory.Cache()
    gmFactory.construct(gmFactoryCache, gContext, gmTree, gmFactoryOptions)

    # Output the surface to an OBJ file
    segments = 720
    if args.output_obj:
        ssurfaces = [ ss[1] for ss in gmFactoryCache.sensitiveSurfaces ]
        acts.examples.writeSurfacesObj(ssurfaces, gContext, [75, 220, 100], segments, "geomodel.obj")

    return

if "__main__" == __name__:
    main()
