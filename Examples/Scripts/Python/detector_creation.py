#!/usr/bin/env python3
import os
import acts
from acts import examples, logging, GeometryContext
from acts.examples.dd4hep import (
    DD4hepDetector,
    DD4hepDetectorOptions,
    DD4hepGeometryService,
)

import json

from common import getOpenDataDetectorDirectory


if "__main__" == __name__:
    odd_xml = getOpenDataDetectorDirectory() / "xml" / "OpenDataDetector.xml"

    print("Using the following xml file: ", odd_xml)

    # Create the dd4hep geometry service and detector
    dd4hepConfig = DD4hepGeometryService.Config()
    dd4hepConfig.logLevel = acts.logging.INFO
    dd4hepConfig.xmlFileNames = [str(odd_xml)]
    dd4hepGeometryService = DD4hepGeometryService(dd4hepConfig)
    dd4hepDetector = DD4hepDetector(dd4hepGeometryService)

    cOptions = DD4hepDetectorOptions(logLevel=acts.logging.VERBOSE, emulateToGraph="")

    # Uncomment if you want to use the geometry id mapping
    # This map can be produced with the 'geometry.py' script
    geoIdMappingFile = None  # 'odd-dd4hep-geoid-mapping.json'
    if geoIdMappingFile is not None:
        # Load the geometry id mapping json file
        with open(geoIdMappingFile) as f:
            # load the file as is
            geometry_id_mapping = json.load(f)
            # create a dictionary with GeometryIdentifier as value
            geometry_id_mapping_patched = {
                int(k): acts.GeometryIdentifier(int(v))
                for k, v in geometry_id_mapping.items()
            }
            # patch the options struct
            acts.examples.dd4hep.attachDD4hepGeoIdMapper(
                cOptions, geometry_id_mapping_patched
            )

    # Context and options
    geoContext = acts.GeometryContext()
    [detector, contextors, store] = dd4hepDetector.finalize(geoContext, cOptions)

    acts.examples.writeDetectorToJsonDetray(geoContext, detector, "odd-detray.json")
