######################################################
# AtlasExtrapolationEngine module
#
# it inherits from Trk__Extrapolator and uses
# the AtlasTrackingGeometrySvc
#
######################################################

# import the include statement
from AthenaCommon.Include import Include, IncludeError, include

# import the ExtrapolationEngine configurable
from ExtrapolationEngine.ExtrapolationEngineConf import Trk__ExtrapolationEngine as ExEngine

# define the class
class AtlasExtrapolationEngine( ExEngine ):
    # constructor
    def __init__(self,name = 'Extrapolation', nameprefix = 'Atlas', ToolOutputLevel = None, TrackingGeometrySvc = None):
       
        # get the correct TrackingGeometry setup
        if not TrackingGeometrySvc :
            from GeometryServices.AtlasTrackingGeometrySvc import AtlasTrackingGeometrySvc
            from AthenaCommon.AppMgr import ServiceMgr as svcMgr
            AtlasTrackingGeometrySvc = svcMgr.AtlasTrackingGeometrySvc 
        else :
            AtlasTrackingGeometrySvc = TrackingGeometrySvc

        # import the ToolSvc
        from AthenaCommon.AppMgr import ToolSvc
        if 'ToolSvc' not in dir() :
            ToolSvc = ToolSvc()

        # load the RungeKutta Propagator
        from TrkExRungeKuttaPropagator.TrkExRungeKuttaPropagatorConf import Trk__RungeKuttaPropagator
        RungeKuttaPropagator = Trk__RungeKuttaPropagator(name = nameprefix+'RungeKuttaPropagator')
        if ToolOutputLevel :
            RungeKuttaPropagator.OutputLevel     = ToolOutputLevel
        ToolSvc += RungeKuttaPropagator
        
        # from the Propagator create a Propagation engine to handle path length
        from ExtrapolationEngine.ExtrapolationEngineConf import Trk__PropagationEngine
        StaticPropagator = Trk__PropagationEngine(name = nameprefix+'StaticPropagation')
        # give the tools it needs 
        StaticPropagator.Propagator               = RungeKuttaPropagator
        # configure output formatting               
        StaticPropagator.OutputPrefix             = '[SP] - '
        StaticPropagator.OutputPostfix            = ' - '
        if ToolOutputLevel : 
            StaticPropagator.OutputLevel          = ToolOutputLevel
        # add to tool service
        ToolSvc += StaticPropagator
        
        # load the material effects engine
        from ExtrapolationEngine.ExtrapolationEngineConf import Trk__MaterialEffectsEngine
        MaterialEffectsEngine = Trk__MaterialEffectsEngine(name = nameprefix+'MaterialEffects')
        # configure output formatting               
        MaterialEffectsEngine.OutputPrefix        = '[ME] - '
        MaterialEffectsEngine.OutputPostfix       = ' - '
        if ToolOutputLevel : 
            MaterialEffectsEngine.OutputLevel     = ToolOutputLevel
        # add to tool service
        ToolSvc += MaterialEffectsEngine
        
        # load the static navigation engine
        from ExtrapolationEngine.ExtrapolationEngineConf import Trk__StaticNavigationEngine
        StaticNavigator = Trk__StaticNavigationEngine(name = nameprefix+'StaticNavigation')
        # give the tools it needs 
        StaticNavigator.PropagationEngine        = StaticPropagator
        StaticNavigator.MaterialEffectsEngine    = MaterialEffectsEngine
        # Geometry name
        StaticNavigator.TrackingGeometry         = AtlasTrackingGeometrySvc.TrackingGeometryName
        # configure output formatting               
        StaticNavigator.OutputPrefix             = '[SN] - '
        StaticNavigator.OutputPostfix            = ' - '
        if ToolOutputLevel : 
            StaticNavigator.OutputLevel              = ToolOutputLevel
        # add to tool service
        ToolSvc += StaticNavigator
        
        
        # load the Static ExtrapolationEngine
        from ExtrapolationEngine.ExtrapolationEngineConf import Trk__StaticEngine
        StaticExtrapolator = Trk__StaticEngine(name = nameprefix+'StaticExtrapolation')
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
        ToolSvc += StaticExtrapolator
       
        # call the base class constructor
        ExEngine.__init__(self, name=nameprefix+'Extrapolation',\
                          ExtrapolationEngines   = [ StaticExtrapolator ], \
                          PropagationEngine      = StaticPropagator, \
                          NavigationEngine       = StaticNavigator, \
                          TrackingGeometrySvc    = AtlasTrackingGeometrySvc, \
                          OutputPrefix           = '[ME] - ', \
                          OutputPostfix          = ' - ')
        # set the output level
        if ToolOutputLevel :
            self.OutputLevel = ToolOutputLevel
