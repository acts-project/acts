from acts import logging, GeometryContext
import json


def setupArgParser():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "--geometryJSON",
        type=str,
        help="JSON dump of the Acts Tracking Geometry",
        required=True,
    )

    return parser


if __name__ == "__main__":
    args = setupArgParser().parse_args()

    # u = acts.UnitConstants
    tgContext = GeometryContext.dangerouslyDefaultConstruct()
    logLevel = logging.INFO

    trackingGeometry = None

    with open(args.geometryJSON, "r") as f:

        from acts.json import TrackingGeometryJsonConverter

        converter = TrackingGeometryJsonConverter()
        geo_json = json.load(f)

        rebuilt_geometry = converter.fromJson(tgContext, args.geometryJSON)
