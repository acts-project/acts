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
        SupportPropertiesB  = MaterialProperties(SupportMaterial, 1., -1.)
        SupportPropertiesC  = MaterialProperties(SupportMaterial, 1., 0.)
        
        
        # the pixel modules
        PixelModuleSmall = DetectorModule(None,8.4,32.0,0.125, Silicon)
        PixelModuleBig   = DetectorModule(None,12.2,32.0,0.15, Silicon)
                
        # the first layer
        PixelLayer0 = CylinderLayer(PixelModuleSmall, 29., 15, 13, 0.2, 2., 0.5, 2., SupportPropertiesA)
        PixelLayer1 = CylinderLayer(PixelModuleSmall, 55., 24, 13, 0.2, 2., 0.5, 2., SupportPropertiesA)
        PixelLayer2 = CylinderLayer(PixelModuleSmall, 88., 40, 13, 0.2, 2., 0.5, 2., SupportPropertiesA)
        PixelLayer3 = CylinderLayer(PixelModuleSmall, 120., 62, 13, 0.2, 2., 0.5, 2., SupportPropertiesB)
        PixelLayer4 = CylinderLayer(PixelModuleBig, 160., 48, 13, 0.2, 2., 0.5, 2., SupportPropertiesB)
        PixelLayer5 = CylinderLayer(PixelModuleBig, 200., 68, 13, 0.2, 2., 0.5, 2., SupportPropertiesB)
        
        # define the pixel barrel volume
        PixelBarrel = BarrelVolume( [ PixelLayer0, PixelLayer1, PixelLayer2, PixelLayer3, PixelLayer4, PixelLayer5 ] ) 

        # lets build some endcap disks
        PixelRing0   = DiscRing(PixelModuleSmall, 65., 24, 0.5)
        PixelRing1   = DiscRing(PixelModuleSmall, 120., 48, 0.5)
        PixelRing2   = DiscRing(PixelModuleSmall, 180, 78, 0.5)
        
        PixelDisc0  = DiscLayer( [ PixelRing0, PixelRing1, PixelRing2 ], 500., 3., 5., SupportPropertiesA)
        PixelDisc1  = DiscLayer( [ PixelRing0, PixelRing1, PixelRing2 ], 580., 3., 5., SupportPropertiesA)
        PixelDisc2  = DiscLayer( [ PixelRing0, PixelRing1, PixelRing2 ], 650., 3., 5., SupportPropertiesA)
        PixelDisc3  = DiscLayer( [ PixelRing0, PixelRing1, PixelRing2 ], 780., 3., 5., SupportPropertiesA)
        # define the pixel endcap volume
        PixelEndcap = EndcapVolume( [ PixelDisc0, PixelDisc1, PixelDisc2, PixelDisc3 ] )
        
        
        # the strip modules 
        StripModuleSmall = DetectorModule(None,16.4,50.5,0.200, Silicon)
        StrupModuleLong  = DetectorModule(None,16.4,50.5,0.200, Silicon)
        
        StripLayer0 = CylinderLayer(StripModuleSmall, 250., 60, 13, -0.2, 5., 5., 2., SupportPropertiesC, -0.02, 0.02, 1.5)
        StripLayer1 = CylinderLayer(StripModuleSmall, 350., 72, 13, -0.2, 5., 5., 2., SupportPropertiesC, -0.02, 0.02, 1.5)
        StripLayer2 = CylinderLayer(StrupModuleLong, 500., 88, 13, -0.2, 5., 5., 2., SupportPropertiesC, -0.02, 0.02, 1.5)
        StripLayer3 = CylinderLayer(StrupModuleLong, 700., 110, 13, -0.2, 5., 5., 2., SupportPropertiesC, -0.02, 0.02, 1.5)
        StripLayer4 = CylinderLayer(StrupModuleLong, 900, 120, 13, -0.2, 5., 5., 2., SupportPropertiesC, -0.02, 0.02, 1.5)

        # define the pixel barrel volume
        StripBarrel = BarrelVolume( [ StripLayer0, StripLayer1, StripLayer2, StripLayer3, StripLayer4 ] ) 


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
        from GenericGeometryTools.GenericGeometryToolsConf import Ats__GenericLayerBuilder as GenericLayerBuilder
        
        
        # # a Pixel layer builder
        PixelLayerBuilder = GenericLayerBuilder('PixelLayerBuilder')
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


        StripLayerBuilder = GenericLayerBuilder('StripLayerBuilder')
        # the ID
        StripLayerBuilder.LayerIdentification               = 'Strip'
        # define the pixel barrel                            
        StripLayerBuilder.CentralLayerRadii                 = StripBarrel.layerRadii()
        StripLayerBuilder.CentralLayerEnvelopeZ             = StripBarrel.layerEnvelopesZ()        
        StripLayerBuilder.CentralLayerMaterialConcentration = StripBarrel.layerMaterialConcentration()        
        StripLayerBuilder.CentralLayerMaterialProperties    = StripBarrel.layerMaterialProperties()        
        StripLayerBuilder.CentralLayerModulesPositionPhi    = StripBarrel.layerModulesPositionPhi()        
        StripLayerBuilder.CentralLayerMoudlesTiltPhi        = StripBarrel.layerModulesTiltPhi()    
        StripLayerBuilder.CentralLayerModulesPositionZ      = StripBarrel.layerModulesPositionZ() 
        StripLayerBuilder.CentralLayerModuleStaggerZ        = StripBarrel.layerModulesStaggerZ()   
        StripLayerBuilder.CentralLayerModulesHalfX          = StripBarrel.layerModulesHalfX()       
        StripLayerBuilder.CentralLayerModulesHalfY          = StripBarrel.layerModulesHalfY()      
        StripLayerBuilder.CentralLayerModulesThickness      = StripBarrel.layerModulesThickness()  
        StripLayerBuilder.CentralLayerModulesMaterial       = StripBarrel.layerModulesMaterial()  
         
        
        # Output steering
        StripLayerBuilder.OutputLevel = outputLevel
        # pixel layer builder is defined
        ToolSvc += StripLayerBuilder

        # Build the Strip Volume
        StripVolumeBuilder = VolumeBuilder("StripVolume")
        # set the volume name
        StripVolumeBuilder.VolumeName                     = 'Strip'
        # build the volume
        StripVolumeBuilder.CylinderVolumeHelper           =  CylinderVolumeHelper   
        StripVolumeBuilder.LayerBuilder                   =  StripLayerBuilder           
        StripVolumeBuilder.LayerArrayCreator              =  LayerArrayCreator        
        StripVolumeBuilder.LayerEnvelopeR                 =  1.    
        StripVolumeBuilder.LayerEnvelopeZ                 =  10.    
        StripVolumeBuilder.VolumeToBeamPipe               =  False  
        ToolSvc += StripVolumeBuilder                                                    
                                                            

        # Build the TrackingGeometry
        GenericGeometryBuilder = GeometryBuilder('GenericGeometry')
        GenericGeometryBuilder.BeamPipeBuilder        = BeamPipeVolumeBuilder
        GenericGeometryBuilder.TrackingVolumeBuilders = [ PixelVolumeBuilder, StripVolumeBuilder  ]
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
        
