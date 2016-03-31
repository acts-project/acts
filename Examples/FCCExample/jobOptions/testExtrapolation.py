###############################################################
#
# Job options
#
#==============================================================

from Gaudi.Configuration import *
from Configurables import ApplicationMgr, THistSvc

#--------------------------------------------------------------
# Private Application Configuration options
#--------------------------------------------------------------

ExToolOutputLevel       = VERBOSE #VERBOSE # INFO #
ExAlgorithmOutputLevel  = VERBOSE #VERBOSE #

from FCCService.FCCServiceConf import Acts__DD4hepGeometrySvc
DD4hepGeometrySvc            = Acts__DD4hepGeometrySvc("DD4hepGeometrySvc", OutputLevel = VERBOSE)
DD4hepGeometrySvc.Detector   = ["file:Examples/FCCExample/FCCDetector/compact/FCCTracker.xml"]

from DD4hepGeometryTools.DD4hepGeometryToolsConf import Acts__DD4hepCylinderGeometryBuilder
DD4hepGeometryBuilder                   = Acts__DD4hepCylinderGeometryBuilder('DD4hepGeometryBuilder')
DD4hepGeometryBuilder.DD4hepGeometrySvc = DD4hepGeometrySvc
DD4hepGeometryBuilder.OutputLevel       = ExToolOutputLevel

from GeometryServices.GeometryServicesConf import Acts__TrackingGeometrySvc
TrackingGeometrySvc = Acts__TrackingGeometrySvc('TrackingGeometrySvc',GeometryBuilder='Acts::DD4hepCylinderGeometryBuilder')
TrackingGeometrySvc.GeometryBuilder.DD4hepGeometrySvc   = DD4hepGeometrySvc
TrackingGeometrySvc.GeometryBuilder                     = DD4hepGeometryBuilder
TrackingGeometrySvc.TrackingGeometryName                = "FCCGeometry"
TrackingGeometrySvc.GeometryProcessors                  = []
TrackingGeometrySvc.OutputLevel                         = ExToolOutputLevel

# configure the json dumper
from JsonWriters.JsonWritersConf import Acts__ParametersJsonWriter as ParametersWriter
JsonParmatersWriter = ParametersWriter('JsonParmatersWriter')

from GenericExtrapolationEngine import GenExEngine
GenericExtrapolationEngine = GenExEngine(name='ExtrapolationEngine', nameprefix='Generic', ToolOutputLevel=ExToolOutputLevel, TrackingGeometrySvc=TrackingGeometrySvc)
ExtrapolationEngine = GenericExtrapolationEngine.extrapolationEngine()
ExtrapolationEngine.OutputLevel = ExToolOutputLevel

#--------------------------------------------------------------
# Algorithm setup
#--------------------------------------------------------------

# Add top algorithms to be run
from ExtrapolationTest.ExtrapolationTestConf import Acts__ExtrapolationEngineTest
ExtrapolationEngineTest = Acts__ExtrapolationEngineTest('ExtrapolationEngineTest')
# how many tests you want per event
ExtrapolationEngineTest.NumberOfTestsPerEvent   = 100000
# parameters mode: 0 - neutral tracks, 1 - charged particles
ExtrapolationEngineTest.ParametersMode          = 1
# do the full test backwards as well
ExtrapolationEngineTest.BackExtrapolation       = False
# Smear the production vertex - standard primary vertex paramters
ExtrapolationEngineTest.SmearOrigin             = False
ExtrapolationEngineTest.SimgaOriginD0           = 0.015
ExtrapolationEngineTest.SimgaOriginZ0           = 55.6
# pT range for testing
ExtrapolationEngineTest.PtMin                   = 1000
ExtrapolationEngineTest.PtMax                   = 10000
# The test range in Eta
ExtrapolationEngineTest.EtaMin                  = -3.5
ExtrapolationEngineTest.EtaMax                  =  3.5
# Configure how you wanna run
ExtrapolationEngineTest.CollectSensitive        = True
ExtrapolationEngineTest.CollectPassive          = True
ExtrapolationEngineTest.CollectBoundary         = True
ExtrapolationEngineTest.CollectMaterial         = True
ExtrapolationEngineTest.SensitiveCurvilinear    = False
ExtrapolationEngineTest.RobustSearch            = False
# the path limit to test
ExtrapolationEngineTest.PathLimit               = -1.
# give it the engine
ExtrapolationEngineTest.ExtrapolationEngine     = ExtrapolationEngine
# the json writer
ExtrapolationEngineTest.ParametersProcessor     = JsonParmatersWriter
# output formatting
ExtrapolationEngineTest.OutputLevel             = ExAlgorithmOutputLevel

THistSvc=THistSvc('THistSvc')
THistSvc.Output = ["val DATAFILE='ExtrapolationEngineTest.root' TYPE='ROOT' OPT='RECREATE'"]
THistSvc.OutputLevel = INFO

ExternalServices = [THistSvc, DD4hepGeometrySvc, TrackingGeometrySvc]
ExternalServices += GenericExtrapolationEngine.containedServices()


ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               ExtSvc=ExternalServices,
               TopAlg=[ExtrapolationEngineTest])

#==============================================================
#
# End of job options file
#
###############################################################