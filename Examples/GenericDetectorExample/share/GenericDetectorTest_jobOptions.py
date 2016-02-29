###############################################################
#
# Job options 
#
#==============================================================

#--------------------------------------------------------------
# ATLAS default Application Configuration options
#--------------------------------------------------------------

# Use McEventSelector so we can run with AthenaMP
import AthenaCommon.AtlasUnixGeneratorJob

#--------------------------------------------------------------
# Private Application Configuration options
#--------------------------------------------------------------

# the global detflags
from AthenaCommon.DetFlags import DetFlags
DetFlags.ID_setOff()
DetFlags.Calo_setOff()
DetFlags.Muon_setOff()

# import the GenericDetector
from GenericDetectorV1 import GenericDetectorConstruction
GenericDetector = GenericDetectorConstruction(name='GenericDetector', outputLevel=VERBOSE)

# Full job is a list of algorithms
from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from AthenaCommon.AppMgr import ToolSvc

# configure the json dumper
from JsonWriters.JsonWritersConf import Ats__GeometryJsonWriter
JsonDumper = Ats__GeometryJsonWriter('GenericGeometrySonDumper')
ToolSvc += JsonDumper

# add the json dumper
TrackingGeometrySvc = GenericDetector.trackingGeometrySvc()
TrackingGeometrySvc.GeometryProcessors = [ JsonDumper ]
 
# Run the GeometryBuildingTestTest
from GeometryBuildingTest.GeometryBuildingTestConf import Ats__TrackingGeometryTest
TrackingGeometryTest = Ats__TrackingGeometryTest('GenericDetectorTest')
TrackingGeometryTest.TrackingGeometrySvc = TrackingGeometrySvc
job += TrackingGeometryTest

# Number of events to be processed (default is until the end of
# input, or -1, however, since we have no input, a limit needs
# to be set explicitly, here, choose 10)
theApp.EvtMax = 1

from AthenaCommon.AppMgr import ServiceMgr
# output level
ServiceMgr.MessageSvc.OutputLevel  = INFO
# increase the number of letter reserved to the alg/tool name from 18 to 30
ServiceMgr.MessageSvc.Format       = "% F%50W%S%7W%R%T %0W%M"
# to change the default limit on number of message
ServiceMgr.MessageSvc.defaultLimit = 9999999  # all messages


#==============================================================
#
# End of job options file
#
###############################################################
