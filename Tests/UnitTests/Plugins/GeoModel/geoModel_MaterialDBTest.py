#!/usr/bin/env python3

"""
CI script which has just the purpose to ensure that the example GeoModel material database can be parsed 
"""


from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--materialDB", help="Path of the material DB to pass", required = True)

args = parser.parse_args()

from acts.geomodel import loadGeoModelMaterialDBfromJSON

if not loadGeoModelMaterialDBfromJSON(args.materialDB):
    print (f"Parsing of the material DB failed {args.materialDB}")
    exit(1)

