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

    p.add_argument("-o", "--output", type=str, default="GeoModel", help="Output file(s) base name")

    p.add_argument(
        "-q",
        "--queries",
        type=str,
        nargs="+",
        default="GeoModelXML",
        help="List of Queries for Published full phys volumes",
    )

    p.add_argument(
        "--table-name",
        type=str,
        default="Blueprint",
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
        "--output-svg",
        help="Write the surfaces to SVG files",
        action="store_true",
        default=False,
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
    ]
    gmFactory = gm.GeoModelDetectorSurfaceFactory(gmFactoryConfig, logging.VERBOSE)
    # The options
    gmFactoryOptions = gm.GeoModelDetectorSurfaceFactory.Options()
    gmFactoryOptions.queries = args.queries
    # The Cache & construct call
    gmFactoryCache = gm.GeoModelDetectorSurfaceFactory.Cache()
    gmFactory.construct(gmFactoryCache, gContext, gmTree, gmFactoryOptions)

    # All surfaces from GeoModel
    gmSurfaces = [ ss[1] for ss in gmFactoryCache.sensitiveSurfaces ]

    # Construct the building hierarchy
    gmBlueprintConfig = gm.GeoModelBlueprintCreater.Config()
    gmBlueprintConfig.detectorSurfaces = gmSurfaces
    gmBlueprintConfig.kdtBinning = [acts.Binning.z, acts.Binning.r]

    gmBlueprintOptions = gm.GeoModelBlueprintCreater.Options()
    gmBlueprintOptions.table = args.table_name
    gmBlueprintOptions.topEntry = args.top_node
    if len(args.top_node_bounds) > 0:
        gmBlueprintOptions.topBoundsOverride = args.top_node_bounds

    gmBlueprintCreater = gm.GeoModelBlueprintCreater(gmBlueprintConfig, logging.VERBOSE)
    gmBlueprint = gmBlueprintCreater.create(gContext, gmTree, gmBlueprintOptions)

    gmCylindricalBuilder = gmBlueprint.convertToBuilder(logging.VERBOSE)

    # Top level geo id generator
    gmGeoIdConfig = GeometryIdGenerator.Config()
    gmGeoIdGenerator = GeometryIdGenerator(
        gmGeoIdConfig, "GeoModelGeoIdGenerator", logging.INFO
    )

    # Create the detector builder
    gmDetectorConfig = DetectorBuilder.Config()
    gmDetectorConfig.name = args.top_node + "_DetectorBuilder"
    gmDetectorConfig.builder = gmCylindricalBuilder
    gmDetectorConfig.geoIdGenerator = gmGeoIdGenerator
    gmDetectorConfig.auxiliary = (
        "GeoModel based Acts::Detector from '" + args.input + "'"
    )

    gmDetectorBuilder = DetectorBuilder(gmDetectorConfig, args.top_node, logging.VERBOSE)
    detector = gmDetectorBuilder.construct(gContext)

    # Output the detector to SVG
    if args.output_svg:
        surfaceStyle = acts.svg.Style()
        surfaceStyle.fillColor = [5, 150, 245]
        surfaceStyle.fillOpacity = 0.5

        surfaceOptions = acts.svg.SurfaceOptions()
        surfaceOptions.style = surfaceStyle

        viewRange = acts.Extent([])
        volumeOptions = acts.svg.DetectorVolumeOptions()
        volumeOptions.surfaceOptions = surfaceOptions

        xyRange = acts.Extent([[acts.Binning.z, [-50, 50]]])
        zrRange = acts.Extent([[acts.Binning.phi, [-0.5, 0.5]]])

        acts.svg.viewDetector(
            gContext,
            detector,
            args.top_node,
            [[ivol, volumeOptions] for ivol in range(detector.numberVolumes())],
            [["xy", ["sensitives", "portals"], xyRange], ["zr", [ "sensitives", "portals", "materials"], zrRange]],
            args.output + "_detector",
        )


    # Output the surface to an OBJ file
    if args.output_obj:
        segments = 720
        ssurfaces = [ss[1] for ss in gmFactoryCache.sensitiveSurfaces]
        acts.examples.writeSurfacesObj(
            ssurfaces, gContext, [75, 220, 100], segments, args.output+"_sensitives.obj"
        )

    return


if "__main__" == __name__:
    main()
