###############################################################
#
# Job options
#
#==============================================================

from Gaudi.Configuration import *
from Configurables import ApplicationMgr

# import the GenericDetector
from GenericDetectorExample.GeometryConstruction import GeometryConstructionGaudi
DetectorGaudi = GeometryConstructionGaudi(name='DetectorGaudi', outputLevel=VERBOSE)

# Run the GeometryBuildingTestTest
from GeometryBuildingTest.GeometryBuildingTestConf import Ats__TrackingGeometryTest
TrackingGeometryTest = Ats__TrackingGeometryTest("TrackingGeometryTest")
TrackingGeometryTest.TrackingGeometrySvc = DetectorGaudi.trackingGeometrySvc()

ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               ExtSvc=[DetectorGaudi.trackingGeometrySvc()],
               TopAlg=[TrackingGeometryTest])


#==============================================================
#
# End of job options file
#
###############################################################

