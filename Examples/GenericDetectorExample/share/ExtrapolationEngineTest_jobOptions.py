##############################################################
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


#--------------------------------------------------------------
# Geometry section
#--------------------------------------------------------------

from AthenaCommon.DetFlags import DetFlags
DetFlags.ID_setOff()
DetFlags.TRT_setOff()
DetFlags.Calo_setOff()
DetFlags.Muon_setOff()

# Full job is a list of algorithms
from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from AthenaCommon.AppMgr import ServiceMgr as svcMgr

from IOVDbSvc.CondDB import conddb
conddb.setGlobalTag('OFLCOND-SIM-00-00-00')

# import the GenericDetector
from GenericDetectorV1 import GenericDetectorConstruction
GenericDetector = GenericDetectorConstruction(name='GenericDetector', outputLevel=VERBOSE)

#--------------------------------------------------------------
# Event related parameters
#--------------------------------------------------------------

# Number of events to be processed (default is until the end of
# input, or -1, however, since we have no input, a limit needs
# to be set explicitly, here, choose 10)
theApp.EvtMax           = 1
ExToolOutputLevel       = INFO #VERBOSE # INFO #
ExAlgorithmOutputLevel  = INFO #VERBOSE #

from AthenaCommon.AppMgr import ServiceMgr
# output level
ServiceMgr.MessageSvc.OutputLevel  = INFO
# increase the number of letter reserved to the alg/tool name from 18 to 30
ServiceMgr.MessageSvc.Format       = "% F%50W%S%7W%R%T %0W%M"
# to change the default limit on number of message
ServiceMgr.MessageSvc.defaultLimit = 9999999  # all messages

#--------------------------------------------------------------
# Tool setup
#--------------------------------------------------------------

# the magnetic field
from MagFieldServices import SetupField
from IOVDbSvc.CondDB import conddb
conddb.addOverride('/GLOBAL/BField/Map','BFieldMap-FullAsym-09-solTil3')

from JsonWriters.JsonWritersConf import Ats__ParametersJsonWriter as ParametersWriter
JsonParmatersWriter = ParametersWriter('JsonParmatersWriter')
ToolSvc += JsonParmatersWriter

from GenericExtrapolationEngine import GenericExtrapolationEngine
ExtrapolationEninge = GenericExtrapolationEngine(name='Extrapolation', nameprefix='Generic', ToolOutputLevel=ExToolOutputLevel, TrackingGeometrySvc=GenericDetector.trackingGeometrySvc())
svcMgr += ExtrapolationEninge

#--------------------------------------------------------------
# Algorithm setup
#--------------------------------------------------------------

# Add top algorithms to be run
from ExtrapolationTest.ExtrapolationTestConf import Ats__ExtrapolationEngineTest
ExtrapolationEngineTest = Ats__ExtrapolationEngineTest('ExtrapolationEngineTest')
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
ExtrapolationEngineTest.PtMin                   = 100
ExtrapolationEngineTest.PtMax                   = 10000
# The test range in Eta                      
ExtrapolationEngineTest.EtaMin                  = -0.05
ExtrapolationEngineTest.EtaMax                  =  0.05
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
ExtrapolationEngineTest.ExtrapolationEngine     = ExtrapolationEninge
# the json writer
ExtrapolationEngineTest.ParametersProcessor     = JsonParmatersWriter
# output formatting
ExtrapolationEngineTest.OutputLevel             = ExAlgorithmOutputLevel
job += ExtrapolationEngineTest   # 1 alg, named 'ExtrapolationEngineTest'


#################################################################
theApp.Dlls += [ 'RootHistCnv' ]
theApp.HistogramPersistency = 'ROOT'

# --- load AuditorSvc
from AthenaCommon.ConfigurableDb import getConfigurable
# --- write out summary of the memory usage
#   | number of events to be skip to detect memory leak
#   | 20 is default. May need to be made larger for complete jobs.
ServiceMgr.AuditorSvc += getConfigurable('ChronoAuditor')()
# --- write out a short message upon entering or leaving each algorithm
#
theApp.AuditAlgorithms = True
theApp.AuditServices   = True
# 
# --- Display detailed size and timing statistics for writing and reading
ServiceMgr.AthenaPoolCnvSvc.UseDetailChronoStat = True


if not hasattr(ServiceMgr, 'THistSvc'):
       from GaudiSvc.GaudiSvcConf import THistSvc
       ServiceMgr += THistSvc()
# add the G4 validation output stream
ServiceMgr.THistSvc.Output += [ "val DATAFILE='ExtrapolationEngineTest.root' TYPE='ROOT' OPT='RECREATE'" ]

#==============================================================
#
# End of job options file
#
###############################################################
