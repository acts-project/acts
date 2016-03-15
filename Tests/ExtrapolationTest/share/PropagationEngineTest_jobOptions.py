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
from AthenaCommon.AppMgr import ToolSvc

from IOVDbSvc.CondDB import conddb
conddb.setGlobalTag('OFLCOND-SIM-00-00-00')


from JsonWriters.JsonWritersConf import Acts__ParametersJsonWriter as ParametersWriter
JsonParmatersWriter = ParametersWriter('JsonParmatersWriter')
ToolSvc += JsonParmatersWriter

#--------------------------------------------------------------
# Event related parameters
#--------------------------------------------------------------

# Number of events to be processed (default is until the end of
# input, or -1, however, since we have no input, a limit needs
# to be set explicitly, here, choose 10)
theApp.EvtMax           = 1

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

from AthenaCommon.CfgGetter import getService
MagFieldSvc =  getService('AtlasFieldSvc')

from RungeKuttaEngine.RungeKuttaEngineConf import Acts__RungeKuttaEngine
RungeKutteEngine = Acts__RungeKuttaEngine(name='RungeKuttaEngine')
RungeKutteEngine.MagneticFieldSvc = MagFieldSvc
ServiceMgr += RungeKutteEngine

#--------------------------------------------------------------
# Algorithm setup
#--------------------------------------------------------------

# Add top algorithms to be run
from ExtrapolationTest.ExtrapolationTestConf import Acts__PropagationEngineTest
PropagationEngineTest = Acts__PropagationEngineTest('PropagationEngineTest')
# give it the engine
PropagationEngineTest.PropagationEngine       = RungeKutteEngine
# the json writer
PropagationEngineTest.ParametersProcessor     = JsonParmatersWriter
# how many tests you want per event 
PropagationEngineTest.NumberOfTestsPerEvent   = 100
# the surface test
PropagationEngineTest.DestinationRadii        = [ 30., 50., 100., 250., 300., 400., 500., 600., 700., 1000. ]
# parameters mode: 0 - neutral tracks, 1 - charged particles 
PropagationEngineTest.ParametersMode          = 1
# do the full test backwards as well            
PropagationEngineTest.EmulatePlaneSurfaces    = True
PropagationEngineTest.ReturnCurvilinear       = False
PropagationEngineTest.BackPropagation         = False
# pT range for testing                        
PropagationEngineTest.PtMin                   = 500
PropagationEngineTest.PtMax                   = 5000
# The test range in Eta                      
PropagationEngineTest.EtaMin                  = -2.5
PropagationEngineTest.EtaMax                  =  2.5
# Configure how you wanna run                  
PropagationEngineTest.PathLimit               = -1.
# output formatting
PropagationEngineTest.OutputLevel             = INFO
job += PropagationEngineTest   # 1 alg, named 'PropagationEngineTest'


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
ServiceMgr.THistSvc.Output += [ "val DATAFILE='PropagationEngineTest.root' TYPE='ROOT' OPT='RECREATE'" ]

#==============================================================
#
# End of job options file
#
###############################################################
