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


def main():
    p = argparse.ArgumentParser()

    p.add_argument("-i", "--input", type=str, default="", help="Input SQL file")

    p.add_argument(
        "-o", "--output", type=str, default="GeoModel", help="Output file(s) base name"
    )

    p.add_argument(
        "-q",
        "--queries",
        type=str,
        nargs="+",
        default="GeoModelXML",
        help="List of Queries for Published full phys volumes",
    )

    p.add_argument(
        "-n",
        "--name-list",
        type=str,
        nargs="+",
        default=[],
        help="List of Name List for the Surface Factory",
    )

    p.add_argument(
        "-ml",
        "--material-list",
        type=str,
        nargs="+",
        default=[],
        help="List of Material List for the Surface Factory",
    )

    p.add_argument(
        "--table-name",
        type=str,
        default="ActsBlueprint",
        help="Name of the blueprint table",
    )

    p.add_argument(
        "-t",
        "--top-node",
        type=str,
        default="",
        help="Name of the top node in the blueprint tree",
    )

    p.add_argument(
        "-b",
        "--top-node-bounds",
        type=str,
        default="",
        help="Table entry string overriding the top node bounds",
    )

    p.add_argument(
        "-m", "--map", type=str, default="", help="Input file for the material map"
    )

    p.add_argument(
        "--convert-subvols",
        help="Convert the children of the top level full phys vol",
        action="store_true",
        default=False,
    )

    p.add_argument(
        "--output-svg",
        help="Write the surfaces to SVG files",
        action="store_true",
        default=False,
    )

    p.add_argument(
        "--output-internals-svg",
        help="Write the internal navigation to SVG files",
        action="store_true",
        default=False,
    )

    p.add_argument(
        "--output-obj",
        help="Write the surfaces to OBJ files",
        action="store_true",
        default=False,
    )

    p.add_argument(
        "--output-json",
        help="Write the surfaces to OBJ files",
        action="store_true",
        default=False,
    )

    p.add_argument(
        "--enable-blueprint",
        help="Enable the usage of the blueprint",
        action="store_true",
        default=False,
    )

    args = p.parse_args()

    gContext = acts.GeometryContext()
    logLevel = logging.INFO

    materialDecorator = None
    if args.map != "":
        print("Loading a material decorator from file:", args.map)
        materialDecorator = acts.IMaterialDecorator.fromFile(args.map)

    # Read the geometry model from the database
    gmTree = acts.geomodel.readFromDb(args.input)
    gmVolFactoryConfig = gm.GeoModelDetectorVolumeFactory.Config()
    gmVolFactoryConfig.materialList = args.material_list
    gmVolFactoryConfig.nameList = args.name_list
    gmVolFactoryConfig.convertSubVolumes = args.convert_subvols
    gmVolFactory = gm.GeoModelDetectorVolumeFactory(gmVolFactoryConfig, logLevel)
    # The options
    gmVolFactoryOptions = gm.GeoModelDetectorVolumeFactory.Options()
    gmVolFactoryOptions.queries = args.queries
    # The Cache & construct call
    gmVolFactoryCache = gm.GeoModelDetectorVolumeFactory.Cache()
    gmVolFactory.construct(gmVolFactoryCache, gContext, gmTree, gmVolFactoryOptions)

    # All surfaces from GeoModel

    # Output the surface to an OBJ file
    if args.output_obj:
        segments = 720
        gmBoxes = gmVolFactoryCache.boundingBoxes
        print(len(gmBoxes))
        acts.examples.writeVolumesObj(
            gmBoxes,
            gContext,
            [75, 220, 100],
            segments,
            args.output + "_vols.obj",
        )

if "__main__" == __name__:
    main()

