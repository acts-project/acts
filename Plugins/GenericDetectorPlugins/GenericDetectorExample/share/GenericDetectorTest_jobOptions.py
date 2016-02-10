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

# Full job is a list of algorithms
from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from AthenaCommon.AppMgr import ToolSvc

# build the beam pipe
from GeometryTools.GeometryToolsConf import Ats__PassiveLayerBuilder as BeamPipeBuilder
BeamPipeBuilder = BeamPipeBuilder('BeamPipeBuilder')
# specify the beam pipe, 0.8 mm of Beryllium here 
BeamPipeBuilder.CentralLayerRadii           = [ 25. ]    
BeamPipeBuilder.CentralLayerHalflengthZ     = [ 3000. ] 
BeamPipeBuilder.CentralLayerThickness       = [ 0.8 ]
BeamPipeBuilder.CentralLayerMaterialX0      = [ 352.8 ]
BeamPipeBuilder.CentralLayerMaterialL0      = [ 407. ]  
BeamPipeBuilder.CentralLayerMaterialA       = [ 9.012 ]
BeamPipeBuilder.CentralLayerMaterialZ       = [ 4. ]
BeamPipeBuilder.CentralLayerMaterialRho     = [ 1.848e-3 ]
ToolSvc += BeamPipeBuilder
#
import math
# 
from GenericGeometryTools.GenericGeometryToolsConf import Agd__GenericLayerBuilder as GenericLayerBuilder
# a Pixel layer builder
PixelLayerBuilder = GenericLayerBuilder('PixelLayerBuilder')
# innermost layer has 
# - 18 staves @ radius 32, 22 staves @ 50, 28 @ 90, 32 @ 130, 38 @ 180   
# - 20 modules along z +/- 500
PixelModuleHalfLengthX = 16
PixelModuleOverlapY    = 3 
PixelModuleHalfLengthY = (500 - 9 * PixelModuleOverlapY)/20.

PixelLayerBuilder.CentralLayerRadii             = [ 32.0, 50., 90., 130., 180. ]
PixelLayerBuilder.CentralLayerEnvelopeZ         = [ 5., 5., 5., 5., 5.]       
PixelLayerBuilder.CentralLayerModulesPhi        = [ ]     
PixelLayerBuilder.CentralLayerMoudlesTiltPhi    = [ 0.2, 0.2, 0.2, 0.2, 0.2 ] 
PixelLayerBuilder.CentralLayerModulesPositionZ  = [ ]
PixelLayerBuilder.CentralLayerModuleStaggerZ    = [ 0.25, 0.25, 0.25, 0.25, 0.25 ]
PixelLayerBuilder.CentralLayerModuleHalfX       = [ PixelModuleHalfLengthX, PixelModuleHalfLengthX, PixelModuleHalfLengthX, PixelModuleHalfLengthX, PixelModuleHalfLengthX ]   
PixelLayerBuilder.CentralLayerModuleHalfY       = [ PixelModuleHalfLengthY, PixelModuleHalfLengthY, PixelModuleHalfLengthY, PixelModuleHalfLengthY, PixelModuleHalfLengthY ] 
PixelLayerBuilder.CentralLayerModuleThickness   = [ 0.125, 0.125, 0.125, 0.125, 0.125 ]



ToolSvc += PixelLayerBuilder





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
