##############################################################
#
# Job options 
#
#==============================================================

#--------------------------------------------------------------
# Gaudi default Application Configuration options
#--------------------------------------------------------------

from Gaudi.Configuration import *
from Configurables import ApplicationMgr, THistSvc

#--------------------------------------------------------------
# Private Application Configuration options
#--------------------------------------------------------------

ExToolOutputLevel       = VERBOSE #VERBOSE # INFO #
ExAlgorithmOutputLevel  = VERBOSE #VERBOSE #
#--------------------------------------------------------------
# Geometry section
#--------------------------------------------------------------

# import the GenericDetector
from GenericDetectorV2Gaudi import GenericDetectorConstructionGaudi
GenericDetector = GenericDetectorConstructionGaudi(name='GenericDetector', outputLevel=VERBOSE)

#--------------------------------------------------------------
# Tool setup
#--------------------------------------------------------------

from JsonWriters.JsonWritersConf import Acts__ParametersJsonWriter as ParametersWriter
JsonParmatersWriter = ParametersWriter('JsonParmatersWriter')

#from ExtrapolationEngine.ExtrapolationEngineConf import Acts__ExtrapolationEngine as ExEngine
#GenericExtrapolationEngine = ExEngine('GenericExtrapolationEngine')

from GenericExtrapolationEngineGaudi import GenericExtrapolationEngineGaudi
GenericExtrapolationEngine = GenericExtrapolationEngineGaudi(name='ExtrapolationEngine', nameprefix='Generic', ToolOutputLevel=ExToolOutputLevel, TrackingGeometrySvc=GenericDetector.trackingGeometrySvc())
ExtrapolationEngine = GenericExtrapolationEngine.extrapolationEngine()

#--------------------------------------------------------------
# Algorithm setup
#--------------------------------------------------------------

# Add top algorithms to be run
from ExtrapolationTest.ExtrapolationTestConf import Acts__ExtrapolationEngineTest
ExtrapolationEngineTest = Acts__ExtrapolationEngineTest('ExtrapolationEngineTest')
# how many tests you want per event 
ExtrapolationEngineTest.NumberOfTestsPerEvent   = 10000
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


#################################################################

THistSvc=THistSvc('THistSvc')
THistSvc.Output = ["val DATAFILE='ExtrapolationEngineTest.root' TYPE='ROOT' OPT='RECREATE'"]
THistSvc.OutputLevel = INFO

ExternalServices = GenericExtrapolationEngine.containedServices()
ExternalServices += [THistSvc]

ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               ExtSvc=ExternalServices,
               TopAlg=[ExtrapolationEngineTest])

#==============================================================
#
# End of job options file
#
###############################################################
