class GeometryConstructionGaudi(object):
    """docstring for GenericDetector"""
    def __init__(self, name, outputLevel):
        super(GeometryConstructionGaudi, self).__init__()
        self.__name__ = name
        self.__tgsvc__ = None
        
        #        from GeometryTools.GeometryToolsConf import Ats__PassiveLayerBuilder as LayerBuilder
        #BeamPipeLayerBuilder = LayerBuilder('BeamPipeLayerBuilder')
        #       from GeometryTools.GeometryToolsConf import Ats__CylinderVolumeHelper as VolumeHelper
        # The Cylinder Volume Helper
        #CylinderVolumeHelper = VolumeHelper('CylinderVolumeHelper')
        
        # Build the TrackingGeometry
        #from GeometryTools.GeometryToolsConf import Ats__CylinderVolumeBuilder as VolumeBuilder
        #BeamPipeVolumeBuilder = VolumeBuilder('BeamPipeVolumeBuilder')
        #BeamPipeVolumeBuilder.CylinderVolumeHelper = CylinderVolumeHelper
        # BeamPipeVolumeBuilder.LayerBuilder = BeamPipeLayerBuilder

        from ToolTest.ToolTestConf import ToolTest2
        Tool2 = ToolTest2("Tool2")

        from ToolTest.ToolTestConf import ToolTest
        Tool1 = ToolTest("Tool1")
        Tool1.ToolTest2 = Tool2

        from GeometryTools.GeometryToolsConf import Ats__CylinderGeometryBuilder as GeometryBuilder
        GenericGeometryBuilder = GeometryBuilder('GenericGeometryBuilder')
        #GenericGeometryBuilder.BeamPipeBuilder = BeamPipeVolumeBuilder
        GenericGeometryBuilder.TestTool = Tool1

        # Establish the TrackingGeometrySvc
        from GeometryServices.GeometryServicesConf import Ats__TrackingGeometrySvc
        GenericTrackingGeometrySvc = Ats__TrackingGeometrySvc('GenericTrackingGeometrySvc')
        GenericTrackingGeometrySvc.GeometryBuilder = GenericGeometryBuilder
        GenericTrackingGeometrySvc.TrackingGeometryName = 'GenericTrackingGeometry'
        GenericTrackingGeometrySvc.GeometryProcessors = []
    
        self.__tgsvc__ = GenericTrackingGeometrySvc
    
    def trackingGeometrySvc(self):
        return self.__tgsvc__
    
    def trackingGeometryName(self):
        return self.__name__
