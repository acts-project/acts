###############################################################
#
# Job options
#
#==============================================================

from Gaudi.Configuration import *
from Configurables import ApplicationMgr

# import the GenericDetector
from GenericDetectorExample.GenericDetectorConstructionGaudi import GenericDetectorConstructionGaudi
GenericDetectorGaudi = GenericDetectorConstructionGaudi(name='GenericDetectorGaudi', outputLevel=VERBOSE)

# Establish the TrackingGeometrySvc
#from GeometryServices.GeometryServicesConf import Ats__TrackingGeometrySvc
#GenericTrackingGeometrySvc = Ats__TrackingGeometrySvc('GenericTrackingGeometrySvc')
#GenericTrackingGeometrySvc.GeometryBuilder = GenericGeometryBuilder
#GenericTrackingGeometrySvc.TrackingGeometryName = 'GenericTrackingGeometry'
#GenericTrackingGeometrySvc.GeometryProcessors = []

# Run the GeometryBuildingTestTest
from GeometryBuildingTest.GeometryBuildingTestConf import Ats__TrackingGeometryTest
TrackingGeometryTest = Ats__TrackingGeometryTest("TrackingGeometryTest")
#TrackingGeometryTest.TrackingGeometrySvc = GenericTrackingGeometrySvc
TrackingGeometryTest.TrackingGeometrySvc = GenericDetectorGaudi.trackingGeometrySvc()

ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               ExtSvc=[GenericDetectorGaudi.trackingGeometrySvc()],
               TopAlg=[TrackingGeometryTest])


#==============================================================
#
# End of job options file
#
###############################################################

