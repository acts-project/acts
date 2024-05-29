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
        "--table-name", type=str, default= "Blueprint", help="Name of the blueprint table"
    )

    p.add_argument(
        "-t", "--top-node", type=str, default= "", help="Name of the top node in the blueprint tree"
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
    gmFactoryConfig.shapeConverters = [ gm.GeoBoxConverter(), gm.GeoTrdConverter(), gm.GeoIntersectionAnnulusConverter() ]
    gmFactory = gm.GeoModelDetectorSurfaceFactory(gmFactoryConfig, logging.VERBOSE)
    # The options
    gmFactoryOptions = gm.GeoModelDetectorSurfaceFactory.Options()
    gmFactoryOptions.queries = args.queries
    # The Cache & construct call
    gmFactoryCache = gm.GeoModelDetectorSurfaceFactory.Cache()
    gmFactory.construct(gmFactoryCache, gContext, gmTree, gmFactoryOptions)

    # Construct the building hierarchy
    gmBlueprintConfig = gm.GeoModelBlueprintCreater.Config()
    gmBlueprintOptions = gm.GeoModelBlueprintCreater.Options()
    gmBlueprintOptions.table = args.table_name
    gmBlueprintOptions.topEntry = args.top_node

    gmBlueprintCreater = gm.GeoModelBlueprintCreater(gmBlueprintConfig, logging.VERBOSE)
    gmBlueprintCreater.create(gContext, gmTree, gmBlueprintOptions)

    # Output the surface to an OBJ file
    segments = 720
    if args.output_obj:
        ssurfaces = [ ss[1] for ss in gmFactoryCache.sensitiveSurfaces ]
        acts.examples.writeSurfacesObj(ssurfaces, gContext, [75, 220, 100], segments, "geomodel.obj")

    return

if "__main__" == __name__:
    main()