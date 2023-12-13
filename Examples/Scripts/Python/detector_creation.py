#!/usr/bin/env python3
import os
import acts
from acts import examples, logging, GeometryContext
from acts.examples.dd4hep import (
    DD4hepDetector,
    DD4hepDetectorOptions,
    DD4hepGeometryService,
)


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

    # Context and options
    geoContext = acts.GeometryContext()
    cOptions = DD4hepDetectorOptions(logLevel=acts.logging.INFO, emulateToGraph="")
    [detector, contextors, store] = dd4hepDetector.finalize(geoContext, cOptions)
