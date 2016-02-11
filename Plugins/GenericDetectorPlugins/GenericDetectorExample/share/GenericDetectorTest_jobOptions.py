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

from AthenaCommon.AppMgr import ToolSvc

# lets build the python detector first -----------------------------------------------
from GenericDetectorDefs import *
# the pixel module
PixelModule = DetectorModule(None,8.4,32.0,0.15)
# the first layer
PixelLayer0 = CylinderLayer(PixelModule, 33., 18, 13, 0.2, 2., 0.5, 5.)
PixelLayer1 = CylinderLayer(PixelModule, 55., 21, 13, 0.2, 2., 0.5, 5.)
PixelLayer2 = CylinderLayer(PixelModule, 88., 24, 13, 0.2, 2., 0.5, 5.)
PixelLayer3 = CylinderLayer(PixelModule, 120., 32, 13, 0.2, 2., 0.5, 5.)
PixelLayer4 = CylinderLayer(PixelModule, 150., 28, 13, 0.2, 2., 0.5, 5.)
# define the pixel barrel volume
PixelBarrel = BarrelVolume( [ PixelLayer0, PixelLayer1, PixelLayer2, PixelLayer3, PixelLayer4 ] ) 

# lets build some endcap disks
PixelRing   = DiscRing(PixelModule, 55., 24, 0.5)
PixelDisc0  = DiscLayer( [ PixelRing ], 500., 0., 5.,)
PixelDisc1  = DiscLayer( [ PixelRing ], 580., 0., 5.,)
PixelDisc2  = DiscLayer( [ PixelRing ], 650., 0., 5.,)
PixelDisc3  = DiscLayer( [ PixelRing ], 700., 0., 5.,)
# define the pixel endcap volume
PixelEndcap = EndcapVolume( [ PixelDisc0, PixelDisc1, PixelDisc2, PixelDisc3 ] )

# -------------------------------------------------------------------------------------
# 
# Builder setup for the layers 
#
# build the beam pipe
from GeometryTools.GeometryToolsConf import Ats__PassiveLayerBuilder as LayerBuilder
from GeometryTools.GeometryToolsConf import Ats__CylinderVolumeBuilder as VolumeBuilder
from GeometryTools.GeometryToolsConf import Ats__CylinderVolumeHelper as VolumeHelper
from GeometryTools.GeometryToolsConf import Ats__CylinderGeometryBuilder as GeometryBuilder
from GeometryTools.GeometryToolsConf import Ats__TrackingVolumeArrayCreator as TrackingVolumeArrayCreator
from GeometryTools.GeometryToolsConf import Ats__LayerArrayCreator as LayerArrayCreator

# Layer Array Creator
LayerArrayCreator = LayerArrayCreator('LayerArrayCreator')
ToolSvc += LayerArrayCreator
# Tracking Volume Array Creator
TrackingVolumeArrayCreator = TrackingVolumeArrayCreator('TrackingVolumeArrayCreator')
ToolSvc += TrackingVolumeArrayCreator

# The Cylinder Volume Helper
CylinderVolumeHelper = VolumeHelper('CylinderVolumeHelper')
CylinderVolumeHelper.LayerArrayCreator = LayerArrayCreator
CylinderVolumeHelper.TrackingVolumeArrayCreator = TrackingVolumeArrayCreator
# done, define it
ToolSvc += CylinderVolumeHelper


BeamPipeBuilder = LayerBuilder('BeamPipeBuilder')
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

BeamPipeVolumeBuilder = VolumeBuilder('BeamPipeVolumeBuilder')
BeamPipeVolumeBuilder.LayerBuilder = BeamPipeBuilder
ToolSvc += BeamPipeVolumeBuilder 

#  
from GenericGeometryTools.GenericGeometryToolsConf import Agd__GenericLayerBuilder as GenericLayerBuilder
# # a Pixel layer builder
PixelLayerBuilder = GenericLayerBuilder('PixelLayerBuilder')
# define the pixel barrel
PixelLayerBuilder.CentralLayerRadii             = [ PixelBarrel.layerRadii()             ]
PixelLayerBuilder.CentralLayerEnvelopeZ         = [ PixelBarrel.layerEnvelopesZ()        ]
PixelLayerBuilder.CentralLayerModulesPhi        = [ PixelBarrel.layerModulesPhi()        ]
PixelLayerBuilder.CentralLayerMoudlesTiltPhi    = [ PixelBarrel.layerModulesTiltPhi()    ]
PixelLayerBuilder.CentralLayerModulesPositionZ  = [ PixelBarrel.layerModulesPositionsZ() ]
PixelLayerBuilder.CentralLayerModuleStaggerZ    = [ PixelBarrel.layerModulesStaggerZ()   ]
PixelLayerBuilder.CentralLayerModuleHalfX       = [ PixelBarrel.layerModulesHalfX()      ] 
PixelLayerBuilder.CentralLayerModuleHalfY       = [ PixelBarrel.layerModulesHalfY()      ]
PixelLayerBuilder.CentralLayerModuleThickness   = [ PixelBarrel.layerModulesThickness()  ]
# define the endcap discs                      
PixelLayerBuilder.PosNegLayerPositionZ          = [ PixelEndcap.layerPositionsZ()        ]
PixelLayerBuilder.PosNegLayerEnvelopeR          = [ PixelEndcap.layerEnvelopesR()        ]
PixelLayerBuilder.PosNegLayerModuleRadii        = [ PixelEndcap.layerModulesRadii()      ]
PixelLayerBuilder.PosNegLayerModuleStaggerR     = [ PixelEndcap.layerModulesStaggerR()   ]
PixelLayerBuilder.PosNegLayerModulesPhi         = [ PixelEndcap.layerModulesPhi()        ]
PixelLayerBuilder.PosNegLayerModulesStaggerPhi  = [ PixelEndcap.layerModulesStaggerPhi() ]
PixelLayerBuilder.PosNegLayerModuleMinHalfX     = [ PixelEndcap.layerModulesMinHalfX()   ]
PixelLayerBuilder.PosNegLayerModuleMaxHalfX     = [ PixelEndcap.layerModulesMaxHalfX()   ]
PixelLayerBuilder.PosNegLayerModuleHalfY        = [ PixelEndcap.layerModulesHalfY()      ]
PixelLayerBuilder.PosNegLayerModuleThickness    = [ PixelEndcap.layerModulesThickness()  ]
# pixel layer builder is defined
ToolSvc += PixelLayerBuilder

# Build the Pixel Volume
PixelVolumeBuilder = VolumeBuilder("PixelVome")
PixelVolumeBuilder.CylinderVolumeHelper =   CylinderVolumeHelper   
PixelVolumeBuilder.LayerBuilder         =   PixelLayerBuilder           
PixelVolumeBuilder.LayerArrayCreator    =   LayerArrayCreator        
PixelVolumeBuilder.LayerEnvelope        =   5.*mm      
ToolSvc += PixelVolumeBuilder

# Build the TrackingGeometry
GenericGeometryBuilder = GeometryBuilder('GenericGeometryBuilder')
GenericGeometryBuilder.BeamPipeBuilder        = BeamPipeVolumeBuilder
GenericGeometryBuilder.TrackingVolumeBuilders = [ PixelVolumeBuilder ]
GenericGeometryBuilder.TrackingVolumeHelper   = CylinderVolumeHelper
ToolSvc += GenericGeometryBuilder

# Establish the TrackingGeometrySvc
from GeometryServices.GeometryServicesConf import Ats__TrackingGeometrySvc
GenericTrackingGeometrySvc = Ats__TrackingGeometrySvc('GenericTrackingGeometrySvc')
GenericTrackingGeometrySvc.GeometryBuilder = GenericGeometryBuilder
GenericTrackingGeometrySvc.TrackingGeometryName = 'GenericTrackingGeometry'
GenericTrackingGeometrySvc.GeometryProcessors = []

from AthenaCommon.AppMgr import ServiceMgr as svcMgr
svcMgr += GenericTrackingGeometrySvc

# Full job is a list of algorithms
from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

# Run the GeometryBuildingTestTest
from GeometryBuildingTest.GeometryBuildingTestConf import Ats__TrackingGeometryTest
TrackingGeometryTest = Ats__TrackingGeometryTest('GenericDetectorTest')
TrackingGeometryTest.TrackingGeometrySvc = GenericTrackingGeometrySvc
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
