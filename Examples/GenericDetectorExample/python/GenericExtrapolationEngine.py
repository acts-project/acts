######################################################
# GenericExtrapoliationEngine module
#
# it inherits from Acts__Extrapolator and uses - currently not
# the AtlasTrackingGeometrySvc
#
######################################################


# import the ExtrapolationEngine configurable
#from ExtrapolationEngine.ExtrapolationEngineConf import Acts__ExtrapolationEngine as ExEngine

### Inheritance currently not possible########################################
#class GenericExtrapolationEngine(ExEngine):
##############################################################################
# define the class
class GenExEngine(object):
    # constructor
    def __init__(self,name = 'Extrapolation', nameprefix = 'Atlas', ToolOutputLevel = None, TrackingGeometrySvc = None, Atlas = False):
        super(GenExEngine, self).__init__()
        self.__name__ = name
        if Atlas:
          # AthenaCommon
          from AthenaCommon.Include import Include, IncludeError, include
          from AthenaCommon.AppMgr import ServiceMgr as svcMgr
          from AthenaCommon.AppMgr import ToolSvc
        
          from AthenaCommon.CfgGetter import getService
          MagneticFieldSvc =  getService('AtlasFieldSvc')

        else:
          #MagneticFieldSvc
          from MagneticFieldServices.MagneticFieldServicesConf import Acts__ConstMagneticFieldSvc
          MagneticFieldSvc = Acts__ConstMagneticFieldSvc('MagneticFieldSvc')
          MagneticFieldSvc.MagneticFieldValue = 0.006
        
        # PropagationEngine
        from RungeKuttaEngine.RungeKuttaEngineConf import Acts__RungeKuttaEngine
        StaticPropagator = Acts__RungeKuttaEngine(name = nameprefix+'StaticPropagation')
        StaticPropagator.MagneticFieldSvc         = MagneticFieldSvc        
        
        # configure output formatting               
        StaticPropagator.OutputPrefix             = '[SP] - '
        StaticPropagator.OutputPostfix            = ' - '
        if ToolOutputLevel : 
            StaticPropagator.OutputLevel          = ToolOutputLevel
        if Atlas:
            svcMgr += StaticPropagator
        # load the material effects engine
        from ExtrapolationEngine.ExtrapolationEngineConf import Acts__MaterialEffectsEngine
        MaterialEffectsEngine = Acts__MaterialEffectsEngine(name = nameprefix+'MaterialEffects')
        # configure output formatting               
        MaterialEffectsEngine.OutputPrefix        = '[ME] - '
        MaterialEffectsEngine.OutputPostfix       = ' - '
        if ToolOutputLevel : 
            MaterialEffectsEngine.OutputLevel     = ToolOutputLevel
        if Atlas:
            svcMgr += MaterialEffectsEngine
        # load the static navigation engine
        from ExtrapolationEngine.ExtrapolationEngineConf import Acts__StaticNavigationEngine
        StaticNavigator = Acts__StaticNavigationEngine(name = nameprefix+'StaticNavigation')
        # give the tools it needs 
        StaticNavigator.PropagationEngine        = StaticPropagator
        StaticNavigator.MaterialEffectsEngine    = MaterialEffectsEngine
        # Geometry name
        StaticNavigator.TrackingGeometrySvc       = TrackingGeometrySvc
        # configure output formatting               
        StaticNavigator.OutputPrefix             = '[SN] - '
        StaticNavigator.OutputPostfix            = ' - '
        if ToolOutputLevel : 
            StaticNavigator.OutputLevel          = ToolOutputLevel
        if Atlas:
            svcMgr += StaticNavigator
        # load the Static ExtrapolationEngine
        from ExtrapolationEngine.ExtrapolationEngineConf import Acts__StaticEngine
        StaticExtrapolator = Acts__StaticEngine(name = nameprefix+'StaticExtrapolation')
        # give the tools it needs 
        StaticExtrapolator.PropagationEngine        = StaticPropagator
        StaticExtrapolator.MaterialEffectsEngine    = MaterialEffectsEngine
        StaticExtrapolator.NavigationEngine         = StaticNavigator
        # configure output formatting               
        StaticExtrapolator.OutputPrefix             = '[SE] - '
        StaticExtrapolator.OutputPostfix            = ' - '
        if ToolOutputLevel : 
            StaticExtrapolator.OutputLevel              = ToolOutputLevel
        if Atlas:
            svcMgr += StaticExtrapolator
    ### Inheritance currently not possible########################################
        # call the base class constructor
        #  ExEngine.__init__(self, name=nameprefix+'Extrapolation',\
        #                    ExtrapolationEngines   = [ StaticExtrapolator ], \
        #                    PropagationEngine      = StaticPropagator, \
        #                    NavigationEngine       = StaticNavigator, \
        #                    TrackingGeometrySvc    = TrackingGeometrySvc, \
        #                    OutputPrefix           = '[ME] - ', \
        #                    OutputPostfix          = ' - ')
        # set the output level
        #  if ToolOutputLevel :
        #    self.OutputLevel = ToolOutputLevel

    ##############################################################################
        from ExtrapolationEngine.ExtrapolationEngineConf import Acts__ExtrapolationEngine as ExEngine
        GenericExtrapolationEngine = ExEngine('GenericExtrapolationEngine')
        GenericExtrapolationEngine.TrackingGeometrySvc = TrackingGeometrySvc
        GenericExtrapolationEngine.ExtrapolationEngines   = [ StaticExtrapolator ]
        GenericExtrapolationEngine.PropagationEngine = StaticPropagator
        GenericExtrapolationEngine.NavigationEngine = StaticNavigator
        GenericExtrapolationEngine.OutputPrefix = '[ME] - '
        GenericExtrapolationEngine.OutputPostfix= ' - '
        # set the output level
        if ToolOutputLevel :
            GenericExtrapolationEngine.OutputLevel = ToolOutputLevel
        if Atlas:
            svcMgr+=GenericExtrapolationEngine
        self.__gexen__ = GenericExtrapolationEngine

        self.__services__= [ TrackingGeometrySvc, MagneticFieldSvc, StaticPropagator, MaterialEffectsEngine, StaticNavigator, StaticExtrapolator, GenericExtrapolationEngine ]

    def containedServices(self):
        return self.__services__
    
    def extrapolationEngine(self):
        return self.__gexen__
    ##############################################################################



