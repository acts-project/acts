import acts
import argparse
from acts import logging, GeometryContext
from acts import geomodel as gm

def main():

    p = argparse.ArgumentParser()
    p.add_argument(
        "-i", "--input", type=str, default="", help="Input SQL file"
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
    gmFactory = gm.GeoModelDetectorSurfaceFactory(logging.VERBOSE)
    # The options
    gmFactoryOptions = gm.GeoModelDetectorSurfaceFactory.Options()
    gmFactoryOptions.queries = [ "GeoModelXML" ]
    # The Cache & construct call
    gmFactoryCache = gm.GeoModelDetectorSurfaceFactory.Cache()
    gmFactory.construct(gmFactoryCache, gContext, gmTree, gmFactoryOptions)
    return

if "__main__" == __name__:
    main()