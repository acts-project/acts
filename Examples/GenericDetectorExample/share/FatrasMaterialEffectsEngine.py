######################################################
# FatrasMaterialEffectsEngine module
#
# It needs to configure all the services and tools
#
######################################################

# import the include statement
from AthenaCommon.Include import Include, IncludeError, include
from AthenaCommon.SystemOfUnits import *

# import the FatrasMaterialEffectsEngine configurable
from FatrasServices.FatrasServicesConf import Acts__FatrasMaterialEffectsEngine as MaterialEffectsEngine
# define the class
class FatrasMaterialEffectsEngine( MaterialEffectsEngine ):
    # constructor
    def __init__(self,name = 'MaterialEffects', nameprefix = 'Fatras', ToolOutputLevel = None):
        
        # AthenaCommon
        from AthenaCommon.AppMgr import ServiceMgr as svcMgr
        from AthenaCommon.AppMgr import ToolSvc
                
        from AtlasCorePlugins.AtlasCorePluginsConf import Acts__AtlasRandomNumberSvc
        RndGenSvc = Acts__AtlasRandomNumberSvc("AtlasRandomNumberSvc")
        svcMgr += RndGenSvc
        
        # Energy Loss Sampler
        from FatrasTools.FatrasToolsConf import Acts__EnergyLossSampler
        EnergyLossSampler = Acts__EnergyLossSampler(name = nameprefix+'EnergyLossSampler')
        EnergyLossSampler.RandomNumberService = RndGenSvc
        if ToolOutputLevel : 
            EnergyLossSampler.OutputLevel          = ToolOutputLevel
        ToolSvc+=EnergyLossSampler
        
        # Electron Energy Loss Sampler
        from FatrasTools.FatrasToolsConf import Acts__ElectronEnergyLossSampler
        ElectronEnergyLossSampler = Acts__ElectronEnergyLossSampler(name = nameprefix+'ElectronEnergyLossSampler')
        ElectronEnergyLossSampler.RandomNumberService = RndGenSvc
        ElectronEnergyLossSampler.ScaleFactor = 1.
        if ToolOutputLevel : 
            ElectronEnergyLossSampler.OutputLevel          = ToolOutputLevel
        ToolSvc+=ElectronEnergyLossSampler
        
        # Hadronic Interaction Sampler
        from FatrasTools.FatrasToolsConf import Acts__HadronicInteractionParametricSampler
        HadronicInteractionSampler = Acts__HadronicInteractionParametricSampler(name = nameprefix+'HadronicInteractionSampler')
        HadronicInteractionSampler.RandomNumberService = RndGenSvc
        HadronicInteractionSampler.PhysicsProcessCode = 121
        HadronicInteractionSampler.MinimumHadronicOutEnergy = 50. ## MeV
        HadronicInteractionSampler.ShortenHadIntChain = False
        if ToolOutputLevel : 
            HadronicInteractionSampler.OutputLevel = ToolOutputLevel
        ToolSvc+=HadronicInteractionSampler
        
        from FatrasTools.FatrasToolsConf import Acts__MultipleScatteringSamplerHighland
        MultipleScatteringSampler = Acts__MultipleScatteringSamplerHighland(name = nameprefix+'MultipleScatteringSampler')
        MultipleScatteringSampler.RandomNumberService = RndGenSvc
        MultipleScatteringSampler.MultipleScatteringLogarithmicTermOn = True
        if ToolOutputLevel : 
            MultipleScatteringSampler.OutputLevel = ToolOutputLevel
        ToolSvc+=MultipleScatteringSampler
        
        from FatrasTools.FatrasToolsConf import Acts__PhotonConversionSampler
        PhotonConversionSampler = Acts__PhotonConversionSampler(name = nameprefix+'PhotonConversionSampler')
        PhotonConversionSampler.RandomNumberService = RndGenSvc
        PhotonConversionSampler.PhysicsProcessCode = 14
        PhotonConversionSampler.MinimumChildEnergy = 50. ## MeV
        PhotonConversionSampler.ChildEnergyScaling = 2.
        if ToolOutputLevel : 
            PhotonConversionSampler.OutputLevel = ToolOutputLevel
        ToolSvc+=PhotonConversionSampler
        
        # call the base class constructor
        MaterialEffectsEngine.__init__(self, name=nameprefix+'MaterialEffects',\
                                           OutputPrefix                  = '[FME] - ', \
                                           OutputPostfix                 = ' - ',\
                                           RandomNumberService           = RndGenSvc,\
                                           DoEnergyLoss                  = True, \
                                           EnergyLossSampler             = EnergyLossSampler,\
                                           DoDedicatedElectronEnergyLoss = True, \
                                           ElectronEnergyLossSampler     = ElectronEnergyLossSampler, \
                                           CreateBremPhotons             = True, \
                                           DoMultipleScattering          = True, \
                                           MultipleScatteringSampler     = MultipleScatteringSampler, \
                                           DoConversion                  = True, \
                                           ConversionSampler             = PhotonConversionSampler, \
                                           DoHadronicInteraction         = True, \
                                           HadronicInteractionSampler    = HadronicInteractionSampler, \
                                           DoDecay                       = True, \
                                           DoPositronAnnihilation        = True, \
                                           MomentumCut                   = 50., \
                                           ParametericScattering         = True, \
                                           BremProcessCode               = 3, \
                                           MinimumBremPhotonMomentum     = 50.)
	
	#set the OutputLevel
        if ToolOutputLevel:
            self.OutputLevel = ToolOutputLevel
            
