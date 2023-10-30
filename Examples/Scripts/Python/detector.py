#!/usr/bin/env python3
import os
import acts
import acts.examples
import acts.examples.dd4hep
from common import getOpenDataDetectorDirectory

dd4hepDetector = acts.examples.dd4hep.DD4hepDetector()

odd_path = getOpenDataDetectorDirectory()
odd_xml = odd_path / "xml" / "OpenDataDetector.xml"

gctx = acts.GeometryContext()

detector, detectorStore, decorators = dd4hepDetector.finalize(gctx, [str(odd_xml)], None)

acts.examples.writeDetectorToJsonDetray(gctx, detector, 'odd')
