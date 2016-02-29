class GenericDetectorConstruction(object):
    """docstring for GenericDetector"""
    def __init__(self, name, outputLevel):
        super(GenericDetectorConstruction, self).__init__()
        self.__name__ = name
        self.__tgsvc__ = None
        
        # ToolSvc and ServiceMgr
        from AthenaCommon.AppMgr import ToolSvc
        from AthenaCommon.AppMgr import ServiceMgr as svcMgr
        import AthenaCommon.Logging as log
        
        # lets build the python detector first -----------------------------------------------
        from GenericGeometryTools.GenericDetectorDefs import DetectorModule
        from GenericGeometryTools.GenericDetectorDefs import CylinderLayer
        from GenericGeometryTools.GenericDetectorDefs import DiscRing
        from GenericGeometryTools.GenericDetectorDefs import DiscLayer
        from GenericGeometryTools.GenericDetectorDefs import BarrelVolume
        from GenericGeometryTools.GenericDetectorDefs import EndcapVolume
        
        # the pixel module
        PixelModule = DetectorModule(None,8.4,32.0,0.15)
        # the first layer
        PixelLayer0 = CylinderLayer(PixelModule, 33., 24, 13, 0.2, 2., 0.5, 2.)
        PixelLayer1 = CylinderLayer(PixelModule, 55., 40, 13, 0.2, 2., 0.5, 2.)
        PixelLayer2 = CylinderLayer(PixelModule, 88., 60, 13, 0.2, 2., 0.5, 2.)
        PixelLayer3 = CylinderLayer(PixelModule, 120., 72, 13, 0.2, 2., 0.5, 2.)
        PixelLayer4 = CylinderLayer(PixelModule, 150., 84, 13, 0.2, 2., 0.5, 2.)
        PixelLayer5 = CylinderLayer(PixelModule, 180., 100, 13, 0.2, 2., 0.5, 2.)
        
        # define the pixel barrel volume
        PixelBarrel = BarrelVolume( [ PixelLayer0, PixelLayer1, PixelLayer2, PixelLayer3, PixelLayer4, PixelLayer5 ] ) 

        # lets build some endcap disks
        PixelRing0   = DiscRing(PixelModule, 65., 24, 0.5)
        PixelRing1   = DiscRing(PixelModule, 120., 48, 0.5)
        PixelRing2   = DiscRing(PixelModule, 180, 78, 0.5)
        
        PixelDisc0  = DiscLayer( [ PixelRing0, PixelRing1, PixelRing2 ], 500., 0., 5.,)
        PixelDisc1  = DiscLayer( [ PixelRing0, PixelRing1, PixelRing2 ], 580., 0., 5.,)
        PixelDisc2  = DiscLayer( [ PixelRing0, PixelRing1, PixelRing2 ], 650., 0., 5.,)
        PixelDisc3  = DiscLayer( [ PixelRing0, PixelRing1, PixelRing2 ], 700., 0., 5.,)
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
        LayerArrayCreator.OutputLevel = outputLevel
        ToolSvc += LayerArrayCreator
        # Tracking Volume Array Creator
        TrackingVolumeArrayCreator = TrackingVolumeArrayCreator('TrackingVolumeArrayCreator')
        TrackingVolumeArrayCreator.OutputLevel = outputLevel
        ToolSvc += TrackingVolumeArrayCreator

        # The Cylinder Volume Helper
        CylinderVolumeHelper = VolumeHelper('CylinderVolumeHelper')
        CylinderVolumeHelper.LayerArrayCreator = LayerArrayCreator
        CylinderVolumeHelper.TrackingVolumeArrayCreator = TrackingVolumeArrayCreator
        CylinderVolumeHelper.OutputLevel = outputLevel
        # done, define it
        ToolSvc += CylinderVolumeHelper


        BeamPipeBuilder = LayerBuilder('BeamPipeBuilder')
        # the identification 
        BeamPipeBuilder.LayerIdentification         = 'BeamPipe'
        # specify the beam pipe, 0.8 mm of Beryllium here 
        BeamPipeBuilder.CentralLayerRadii           = [ 20   ]    
        BeamPipeBuilder.CentralLayerHalflengthZ     = [ 200. ] 
        BeamPipeBuilder.CentralLayerThickness       = [ 0.8 ]
        BeamPipeBuilder.CentralLayerMaterialX0      = [ 352.8 ]
        BeamPipeBuilder.CentralLayerMaterialL0      = [ 407. ]  
        BeamPipeBuilder.CentralLayerMaterialA       = [ 9.012 ]
        BeamPipeBuilder.CentralLayerMaterialZ       = [ 4. ]
        BeamPipeBuilder.CentralLayerMaterialRho     = [ 1.848e-3 ]
        # output
        BeamPipeBuilder.OutputLevel = outputLevel
        ToolSvc += BeamPipeBuilder

        BeamPipeVolumeBuilder = VolumeBuilder('BeamPipeVolumeBuilder')
        # set the volume nome
        BeamPipeVolumeBuilder.VolumeName           = BeamPipeBuilder.LayerIdentification 
        # and some parameters
        BeamPipeVolumeBuilder.CylinderVolumeHelper = CylinderVolumeHelper   
        BeamPipeVolumeBuilder.LayerBuilder         = BeamPipeBuilder
        BeamPipeVolumeBuilder.LayerArrayCreator    = LayerArrayCreator        
        BeamPipeVolumeBuilder.LayerEnvelope        = 1.   
        BeamPipeVolumeBuilder.OutputLevel          = outputLevel  
        ToolSvc += BeamPipeVolumeBuilder 
        #  
        from GenericGeometryTools.GenericGeometryToolsConf import Ats__GenericLayerBuilder as GenericLayerBuilder
        # # a Pixel layer builder
        PixelLayerBuilder = GenericLayerBuilder('PixelLayerBuilder')
        # the ID
        PixelLayerBuilder.LayerIdentification             = 'Pixel'
        # define the pixel barrel
        PixelLayerBuilder.CentralLayerRadii               = PixelBarrel.layerRadii()             
        PixelLayerBuilder.CentralLayerEnvelopeZ           = PixelBarrel.layerEnvelopesZ()        
        PixelLayerBuilder.CentralLayerModulesPositionPhi  = PixelBarrel.layerModulesPositionPhi()        
        PixelLayerBuilder.CentralLayerMoudlesTiltPhi      = PixelBarrel.layerModulesTiltPhi()    
        PixelLayerBuilder.CentralLayerModulesPositionZ    = PixelBarrel.layerModulesPositionZ() 
        PixelLayerBuilder.CentralLayerModuleStaggerZ      = PixelBarrel.layerModulesStaggerZ()   
        PixelLayerBuilder.CentralLayerModulesHalfX        = PixelBarrel.layerModulesHalfX()       
        PixelLayerBuilder.CentralLayerModulesHalfY        = PixelBarrel.layerModulesHalfY()      
        PixelLayerBuilder.CentralLayerModulesThickness    = PixelBarrel.layerModulesThickness()  
        # define the endcap discs                      
        PixelLayerBuilder.PosNegLayerPositionZ            = PixelEndcap.layerPositionsZ()       
        PixelLayerBuilder.PosNegLayerEnvelopeR            = PixelEndcap.layerEnvelopesR()       
        PixelLayerBuilder.PosNegLayerModulesRadii         = PixelEndcap.layerModulesRadii()     
        PixelLayerBuilder.PosNegLayerModuleStaggerR       = PixelEndcap.layerModulesStaggerR()  
        PixelLayerBuilder.PosNegLayerModulesInPhi         = PixelEndcap.layerModulesInPhi()       
        PixelLayerBuilder.PosNegLayerModulesPositionPhi   = PixelEndcap.layerModulesPositionPhi()       
        PixelLayerBuilder.PosNegLayerModulesStaggerPhi    = PixelEndcap.layerModulesStaggerPhi()
        PixelLayerBuilder.PosNegLayerModulesMinHalfX      = PixelEndcap.layerModulesMinHalfX() 
        PixelLayerBuilder.PosNegLayerModulesMaxHalfX      = PixelEndcap.layerModulesMaxHalfX() 
        PixelLayerBuilder.PosNegLayerModulesHalfY         = PixelEndcap.layerModulesHalfY()    
        PixelLayerBuilder.PosNegLayerModulesThickness     = PixelEndcap.layerModulesThickness()
        # Output steering
        PixelLayerBuilder.OutputLevel = outputLevel
        # pixel layer builder is defined
        ToolSvc += PixelLayerBuilder

        # Build the Pixel Volume
        PixelVolumeBuilder = VolumeBuilder("PixelVome")
        # set the volume name
        PixelVolumeBuilder.VolumeName                     = 'Pixel'
        # build the volume
        PixelVolumeBuilder.CylinderVolumeHelper           =  CylinderVolumeHelper   
        PixelVolumeBuilder.LayerBuilder                   =  PixelLayerBuilder           
        PixelVolumeBuilder.LayerArrayCreator              =  LayerArrayCreator        
        PixelVolumeBuilder.LayerEnvelope                  =  1.    
        PixelVolumeBuilder.VolumeToBeamPipe               =  False  
        ToolSvc += PixelVolumeBuilder

        # Build the TrackingGeometry
        GenericGeometryBuilder = GeometryBuilder('GenericGeometry')
        GenericGeometryBuilder.BeamPipeBuilder        = BeamPipeVolumeBuilder
        GenericGeometryBuilder.TrackingVolumeBuilders = [ PixelVolumeBuilder ]
        GenericGeometryBuilder.TrackingVolumeHelper   = CylinderVolumeHelper
        GenericGeometryBuilder.OutputLevel            = outputLevel
        ToolSvc += GenericGeometryBuilder

        # Establish the TrackingGeometrySvc
        from GeometryServices.GeometryServicesConf import Ats__TrackingGeometrySvc
        GenericTrackingGeometrySvc = Ats__TrackingGeometrySvc('GenericTrackingGeometrySvc')
        GenericTrackingGeometrySvc.GeometryBuilder = GenericGeometryBuilder
        GenericTrackingGeometrySvc.TrackingGeometryName = 'GenericTrackingGeometry'
        GenericTrackingGeometrySvc.GeometryProcessors = []

        svcMgr += GenericTrackingGeometrySvc
        self.__tgsvc__ = GenericTrackingGeometrySvc
    
    def trackingGeometrySvc(self):
        return self.__tgsvc__
        
    def trackingGeometryName(self):
        return self.__name__
        
