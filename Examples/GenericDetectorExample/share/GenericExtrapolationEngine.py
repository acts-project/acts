######################################################
# GenericExtrapoliationEngine module
#
# it inherits from Ats__Extrapolator and uses
# the AtlasTrackingGeometrySvc
#
######################################################

# import the include statement
from AthenaCommon.Include import Include, IncludeError, include

# import the ExtrapolationEngine configurable
from ExtrapolationEngine.ExtrapolationEngineConf import Ats__ExtrapolationEngine as ExEngine

# define the class
class GenericExtrapolationEngine( ExEngine ):
    # constructor
    def __init__(self,name = 'Extrapolation', nameprefix = 'Atlas', ToolOutputLevel = None, TrackingGeometrySvc = None):
       
        # AthenaCommon
        from AthenaCommon.AppMgr import ServiceMgr as svcMgr
        from AthenaCommon.AppMgr import ToolSvc
        
        from AthenaCommon.CfgGetter import getService
        MagneticFieldSvc =  getService('AtlasFieldSvc')
        
        # PropagationEngine
        from RungeKuttaEngine.RungeKuttaEngineConf import Ats__RungeKuttaEngine
        StaticPropagator = Ats__RungeKuttaEngine(name = nameprefix+'StaticPropagation')
        StaticPropagator.MagneticFieldSvc         = MagneticFieldSvc        
        
        # configure output formatting               
        StaticPropagator.OutputPrefix             = '[SP] - '
        StaticPropagator.OutputPostfix            = ' - '
        if ToolOutputLevel : 
            StaticPropagator.OutputLevel          = ToolOutputLevel
        # add to tool service
        svcMgr += StaticPropagator
        
        # load the material effects engine
        from ExtrapolationEngine.ExtrapolationEngineConf import Ats__MaterialEffectsEngine
        MaterialEffectsEngine = Ats__MaterialEffectsEngine(name = nameprefix+'MaterialEffects')
        # configure output formatting               
        MaterialEffectsEngine.OutputPrefix        = '[ME] - '
        MaterialEffectsEngine.OutputPostfix       = ' - '
        if ToolOutputLevel : 
            MaterialEffectsEngine.OutputLevel     = ToolOutputLevel
        # add to tool service
        svcMgr += MaterialEffectsEngine
        
        # load the static navigation engine
        from ExtrapolationEngine.ExtrapolationEngineConf import Ats__StaticNavigationEngine
        StaticNavigator = Ats__StaticNavigationEngine(name = nameprefix+'StaticNavigation')
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
        # add to tool service
        svcMgr += StaticNavigator
        
        
        # load the Static ExtrapolationEngine
        from ExtrapolationEngine.ExtrapolationEngineConf import Ats__StaticEngine
        StaticExtrapolator = Ats__StaticEngine(name = nameprefix+'StaticExtrapolation')
        # give the tools it needs 
        StaticExtrapolator.PropagationEngine        = StaticPropagator
        StaticExtrapolator.MaterialEffectsEngine    = MaterialEffectsEngine
        StaticExtrapolator.NavigationEngine         = StaticNavigator
        # configure output formatting               
        StaticExtrapolator.OutputPrefix             = '[SE] - '
        StaticExtrapolator.OutputPostfix            = ' - '
        if ToolOutputLevel : 
            StaticExtrapolator.OutputLevel              = ToolOutputLevel
        # add to tool service
        svcMgr += StaticExtrapolator
       
        # call the base class constructor
        ExEngine.__init__(self, name=nameprefix+'Extrapolation',\
                          ExtrapolationEngines   = [ StaticExtrapolator ], \
                          PropagationEngine      = StaticPropagator, \
                          NavigationEngine       = StaticNavigator, \
                          TrackingGeometrySvc    = TrackingGeometrySvc, \
                          OutputPrefix           = '[ME] - ', \
                          OutputPostfix          = ' - ')
        # set the output level
        if ToolOutputLevel :
            self.OutputLevel = ToolOutputLevel
