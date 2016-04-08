###############################################################
#
# Job options
#
#==============================================================

from Gaudi.Configuration import *
from Configurables import ApplicationMgr

from FCCService.FCCServiceConf import Acts__DD4hepGeometrySvc
DD4hepGeometrySvc            = Acts__DD4hepGeometrySvc("DD4hepGeometrySvc", OutputLevel = VERBOSE)
DD4hepGeometrySvc.Detector   = ["file:Examples/FCCExample/FCCDetector/compact/TKLayoutTracker.xml"]

# configure the json dumper
from JsonWriters.JsonWritersConf import Acts__GeometryJsonWriter
JsonDumper = Acts__GeometryJsonWriter('GenericGeometrySonDumper')

from DD4hepGeometryTools.DD4hepGeometryToolsConf import Acts__DD4hepCylinderGeometryBuilder
DD4hepGeometryBuilder                   = Acts__DD4hepCylinderGeometryBuilder('DD4hepGeometryBuilder')
DD4hepGeometryBuilder.DD4hepGeometrySvc = DD4hepGeometrySvc

from GeometryServices.GeometryServicesConf import Acts__TrackingGeometrySvc
TrackingGeometrySvc = Acts__TrackingGeometrySvc('TrackingGeometrySvc',GeometryBuilder='Acts::DD4hepCylinderGeometryBuilder')
TrackingGeometrySvc.GeometryBuilder.DD4hepGeometrySvc   = DD4hepGeometrySvc
TrackingGeometrySvc.GeometryBuilder                     = DD4hepGeometryBuilder
TrackingGeometrySvc.TrackingGeometryName                = "FCCGeometry"
TrackingGeometrySvc.GeometryProcessors                  = [ JsonDumper ]

from GeometryBuildingTest.GeometryBuildingTestConf import Acts__TrackingGeometryTest
TrackingGeometryTest                     = Acts__TrackingGeometryTest('DD4hepDetectorTest')
TrackingGeometryTest.TrackingGeometrySvc = TrackingGeometrySvc


ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               ExtSvc=[DD4hepGeometrySvc, TrackingGeometrySvc],
               TopAlg=[ TrackingGeometryTest])

#==============================================================
#
# End of job options file
#
###############################################################
