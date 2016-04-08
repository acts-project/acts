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

        from GenericGeometryTools.GenericDetectorDefs import Material
        from GenericGeometryTools.GenericDetectorDefs import MaterialProperties
        from GenericGeometryTools.GenericDetectorDefs import DetectorModule
        from GenericGeometryTools.GenericDetectorDefs import CylinderLayer
        from GenericGeometryTools.GenericDetectorDefs import DiscRing
        from GenericGeometryTools.GenericDetectorDefs import DiscLayer
        from GenericGeometryTools.GenericDetectorDefs import BarrelVolume
        from GenericGeometryTools.GenericDetectorDefs import EndcapVolume
        
        # Silicon Material
        Silicon             = Material('Silicon', 95.7, 465.2, 28.03,14.,2.32e-3)
        SupportMaterial     = Material('Support', 95.7, 465.2, 28.03,14.,2.32e-3)
        
        SupportPropertiesA  = MaterialProperties(SupportMaterial, 1.,  1.)
        
        # the pixel modules
        PixelModuleSmall = DetectorModule(None,8.4,32.0,0.15, Silicon)
        PixelModuleBig   = DetectorModule(None,12.2,32.0,0.15, Silicon)
        
        # the first layer
        PixelLayer0 = CylinderLayer(PixelModuleSmall, 29., 15, 13, 0.18, 2., 0.5, 2., SupportPropertiesA)
        PixelLayer1 = CylinderLayer(PixelModuleSmall, 55., 24, 13, 0.18, 2., 0.5, 2., SupportPropertiesA)
        PixelLayer2 = CylinderLayer(PixelModuleSmall, 88., 40, 13, 0.2, 2., 0.5, 2., SupportPropertiesA)
        
        # define the pixel barrel volume
        PixelBarrel = BarrelVolume( [ PixelLayer0, PixelLayer1, PixelLayer2 ] ) 

        # lets build some endcap disks
        PixelRing0   = DiscRing(PixelModuleSmall, 65., 24, 0.5)
        
        PixelDisc0  = DiscLayer( [ PixelRing0 ], 500., 3., 5., SupportPropertiesA)
        PixelDisc1  = DiscLayer( [ PixelRing0 ], 580., 3., 5., SupportPropertiesA)
        PixelDisc2  = DiscLayer( [ PixelRing0 ], 650., 3., 5., SupportPropertiesA)
        # define the pixel endcap volume
        PixelEndcap = EndcapVolume( [ PixelDisc0, PixelDisc1, PixelDisc2 ] )

        # -------------------------------------------------------------------------------------
        # 
        # Builder setup for the layers 
        #
        # import the tools
        from GeometryTools.GeometryToolsConf import Acts__LayerCreator as LayerCreator
        from GeometryTools.GeometryToolsConf import Acts__PassiveLayerBuilder as LayerBuilder
        from GeometryTools.GeometryToolsConf import Acts__CylinderVolumeBuilder as VolumeBuilder
        from GeometryTools.GeometryToolsConf import Acts__CylinderVolumeHelper as VolumeHelper
        from GeometryTools.GeometryToolsConf import Acts__CylinderGeometryBuilder as GeometryBuilder
        from GeometryTools.GeometryToolsConf import Acts__TrackingVolumeArrayCreator as TrackingVolumeArrayCreator
        from GeometryTools.GeometryToolsConf import Acts__LayerArrayCreator as LayerArrayCreator
        from GeometryTools.GeometryToolsConf import Acts__SurfaceArrayCreator as SurfaceArrayCreator

        # Surface Array Creator
        SurfaceArrayCreator = SurfaceArrayCreator('SurfaceArrayCreator')
        SurfaceArrayCreator.OutputLevel = outputLevel
        ToolSvc += SurfaceArrayCreator
        # Layer Creator
        LayerCreator = LayerCreator('LayerCreator')
        LayerCreator.SurfaceArrayCreator = SurfaceArrayCreator
        LayerCreator.OutputLevel = outputLevel
        ToolSvc += LayerCreator
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
        BeamPipeBuilder.CentralLayerRadii           = [ 21. ]    
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
        BeamPipeVolumeBuilder.LayerEnvelopeR       = 1.
        BeamPipeVolumeBuilder.LayerEnvelopeZ       = 1.
        BeamPipeVolumeBuilder.OutputLevel          = outputLevel  
        ToolSvc += BeamPipeVolumeBuilder 
        #  
        from GenericGeometryTools.GenericGeometryToolsConf import Acts__GenericLayerBuilder as GenericLayerBuilder
        # # a Pixel layer builder
        PixelLayerBuilder = GenericLayerBuilder('PixelLayerBuilder')
        # the layer creator
        PixelLayerBuilder.LayerCreator                      = LayerCreator
        # the ID
        PixelLayerBuilder.LayerIdentification               = 'Pixel'
        # define the pixel barrel                            
        PixelLayerBuilder.CentralLayerRadii                 = PixelBarrel.layerRadii()
        PixelLayerBuilder.CentralLayerEnvelopeZ             = PixelBarrel.layerEnvelopesZ()        
        PixelLayerBuilder.CentralLayerMaterialConcentration = PixelBarrel.layerMaterialConcentration()        
        PixelLayerBuilder.CentralLayerMaterialProperties    = PixelBarrel.layerMaterialProperties()        
        PixelLayerBuilder.CentralLayerModulesPositionPhi    = PixelBarrel.layerModulesPositionPhi()        
        PixelLayerBuilder.CentralLayerMoudlesTiltPhi        = PixelBarrel.layerModulesTiltPhi()    
        PixelLayerBuilder.CentralLayerModulesPositionZ      = PixelBarrel.layerModulesPositionZ() 
        PixelLayerBuilder.CentralLayerModuleStaggerZ        = PixelBarrel.layerModulesStaggerZ()   
        PixelLayerBuilder.CentralLayerModulesHalfX          = PixelBarrel.layerModulesHalfX()       
        PixelLayerBuilder.CentralLayerModulesHalfY          = PixelBarrel.layerModulesHalfY()      
        PixelLayerBuilder.CentralLayerModulesThickness      = PixelBarrel.layerModulesThickness()  
        PixelLayerBuilder.CentralLayerModulesMaterial       = PixelBarrel.layerModulesMaterial()  
                                                            
        # define the endcap discs                           
        PixelLayerBuilder.PosNegLayerPositionZ              = PixelEndcap.layerPositionsZ()       
        PixelLayerBuilder.PosNegLayerEnvelopeR              = PixelEndcap.layerEnvelopesR()  
        PixelLayerBuilder.PosNegLayerMaterialConcentration  = PixelEndcap.layerMaterialConcentration()        
        PixelLayerBuilder.PosNegLayerMaterialProperties     = PixelEndcap.layerMaterialProperties()        
        PixelLayerBuilder.PosNegLayerModulesRadii           = PixelEndcap.layerModulesRadii()     
        PixelLayerBuilder.PosNegLayerModuleStaggerR         = PixelEndcap.layerModulesStaggerR()  
        PixelLayerBuilder.PosNegLayerModulesInPhi           = PixelEndcap.layerModulesInPhi()       
        PixelLayerBuilder.PosNegLayerModulesPositionPhi     = PixelEndcap.layerModulesPositionPhi()       
        PixelLayerBuilder.PosNegLayerModulesStaggerPhi      = PixelEndcap.layerModulesStaggerPhi()
        PixelLayerBuilder.PosNegLayerModulesMinHalfX        = PixelEndcap.layerModulesMinHalfX() 
        PixelLayerBuilder.PosNegLayerModulesMaxHalfX        = PixelEndcap.layerModulesMaxHalfX() 
        PixelLayerBuilder.PosNegLayerModulesHalfY           = PixelEndcap.layerModulesHalfY()    
        PixelLayerBuilder.PosNegLayerModulesThickness       = PixelEndcap.layerModulesThickness()
        PixelLayerBuilder.PosNegLayerModulesMaterial        = PixelEndcap.layerModulesMaterial()
        
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
        PixelVolumeBuilder.LayerEnvelopeR                 =  1.    
        PixelVolumeBuilder.LayerEnvelopeZ                 =  10.   
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
        from GeometryServices.GeometryServicesConf import Acts__TrackingGeometrySvc
        GenericTrackingGeometrySvc = Acts__TrackingGeometrySvc('GenericTrackingGeometrySvc')
        GenericTrackingGeometrySvc.GeometryBuilder = GenericGeometryBuilder
        GenericTrackingGeometrySvc.TrackingGeometryName = 'GenericTrackingGeometry'
        GenericTrackingGeometrySvc.GeometryProcessors = []

        svcMgr += GenericTrackingGeometrySvc
        self.__tgsvc__ = GenericTrackingGeometrySvc
    
    def trackingGeometrySvc(self):
        return self.__tgsvc__
        
    def trackingGeometryName(self):
        return self.__name__
        
