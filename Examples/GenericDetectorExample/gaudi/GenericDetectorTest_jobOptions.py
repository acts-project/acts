###############################################################
#
# Job options
#
#==============================================================

from Gaudi.Configuration import *
from Configurables import ApplicationMgr

# import the GenericDetector
from GenericDetectorV2 import GenericDetectorConstruction
GenericDetector = GenericDetectorConstruction(name='GenericDetector', outputLevel=VERBOSE)

# configure the json dumper
from JsonWriters.JsonWritersConf import Acts__GeometryJsonWriter
JsonDumper = Acts__GeometryJsonWriter('GenericGeometrySonDumper')

# add the json dumper
TrackingGeometrySvc = GenericDetector.trackingGeometrySvc()
TrackingGeometrySvc.GeometryProcessors = [ JsonDumper ]

# Run the GeometryBuildingTestTest
from GeometryBuildingTest.GeometryBuildingTestConf import Acts__TrackingGeometryTest
TrackingGeometryTest = Acts__TrackingGeometryTest("TrackingGeometryTest")
TrackingGeometryTest.TrackingGeometrySvc = TrackingGeometrySvc

ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               ExtSvc=[TrackingGeometrySvc],
               TopAlg=[TrackingGeometryTest])

#==============================================================
#
# End of job options file
#
###############################################################

