###############################################################
#
# Job options
#
#==============================================================

from Gaudi.Configuration import *
from Configurables import ApplicationMgr

# import the GenericDetector
from GenericDetectorGaudi import GenericDetectorConstructionGaudi
GenericDetectorGaudi = GenericDetectorConstructionGaudi(name='GenericDetectorGaudi', outputLevel=VERBOSE)
# Run the GeometryBuildingTestTest
from GeometryBuildingTest.GeometryBuildingTestConf import Ats__TrackingGeometryTest
TrackingGeometryTest = Ats__TrackingGeometryTest("TrackingGeometryTest")
TrackingGeometryTest.TrackingGeometrySvc = GenericDetectorGaudi.trackingGeometrySvc()

ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               ExtSvc=[],
               TopAlg=[TrackingGeometryTest])


#==============================================================
#
# End of job options file
#
###############################################################

