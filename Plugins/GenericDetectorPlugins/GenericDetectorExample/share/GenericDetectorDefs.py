import math

""" Detector Module definition """
class DetectorModule(object):
    
    """docstring for DetectorModule: defines a simple detector moduel for generic detetors"""
    def __init__(self, minHalfX, maxHalfX, halfY, thickness):
        super(DetectorModule, self).__init__()
        self.__minHalfX__  = minHalfX
        self.__maxHalfX__  = maxHalfX
        self.__halfY__     = halfY
        self.__thickness__ = thickness

    """return methods"""
    def minHalfX(self):
        return self.__minHalfX__

    """return methods"""
    def maxHalfX(self):
        return self.__maxHalfX__

    """return methods"""
    def halfY(self):
        return self.__halfY__

    """return methods"""
    def thickness(self):
        return self.__thickness__

""" Cylinder Layer definition """
class CylinderLayer(object):
    """docstring for CylinderLayer"""
    def __init__(self, module, radius, numPhi, numZ, tiltPhi, overlapZ, staggerZ, evelopeCover ):
        super(CylinderLayer, self).__init__()
        self.__radius__   =   radius 
        self.__module__   =   module 
        self.__numPhi__   =   numPhi 
        self.__numZ__     =   numZ 
        self.__tiltPhi__  =   tiltPhi 
        self.__overlapZ__ =   overlapZ
        self.__staggerZ__ =   staggerZ
        self.__eCover__   =   evelopeCover
        
    def radius(self):
        return self.__radius__
    
    def module(self):
        return self.__module__
         
    def modulesInPhi(self) :
        return self.__numPhi__
        
    def modulePositionPhi(self) :
        phiPositions = []
        phiStep = 2.*math.pi/self.__numPhi__;
        for iphi in xrange(self.__numPhi__) :
            phiPositions += [ -math.pi+iphi*phiStep]
        return phiPositions
        
    def modulesInZ(self) :    
        return self.__numZ__
        
    def modulePositionZ(self) :
        zPositions  = []
        firstToLast = (self.__numZ__-1)*(self.__module__.halfY()-self.__overlapZ__)
        zStep = firstToLast/(self.__numZ__-1)
        for iz in xrange(self.__numZ__) :
            zPositions += [ -0.5*firstToLast + iz*zStep ]
        return zPositions
        
    def tiltPhi(self) :
        return self.__tiltPhi__
        
    def overlapZ(self) :     
        return self.__overlapZ__
        
    def staggerZ(self) :
        return self.__staggerZ__
        
    def enevlopeCover(self) :
            return self.__eCover__
                        
""" Disc Ring Definition """            
class DiscRing(object):
    """docstring for DiscRing"""
    def __init__(self, module, radius, numPhi, staggerPhi):
        super(DiscRing, self).__init__()
        self.__radius__     =   radius 
        self.__module__     =   module 
        self.__numPhi__     =   numPhi 
        self.__staggerPhi__ =   staggerPhi
        
    def radius(self):
        return self.__radius__
    
    def module(self):
        return self.__module__
         
    def modulesInPhi(self) :
        return self.__numPhi__
        
    def modulePositionPhi(self) :
        phiPositions = []
        phiStep = 2.*math.pi/self.__numPhi__;
        for iphi in xrange(self.__numPhi__) :
            phiPositions += [ -math.pi+iphi*phiStep]
            return phiPositions
            
    def phiStagger(self):
        return self.__staggerPhi__


""" Disc Layer Definition """            
class DiscLayer(object):
    
    """docstring for DiscLayer"""
    def __init__(self, rings, zPosition, rStagger, rEnvelope):
        super(DiscLayer, self).__init__()
        self.__zPosition__       = zPosition
        self.__rStagger__        = rStagger
        self.__rEnvelope__       = rEnvelope
        self.__rings__           = rings
        # the ring components
        self.__rRadii__           = []
        self.__nPhis__            = []
        self.__rphistaggers__     = []
        self.__minHalfLengtshX__  = []
        self.__maxHalfLengtshX__  = []
        self.__halfLengtshX__     = []
        self.__moduleThickness__  = []
        
        for ring in self.__rings__ :
            self.__rRadii__          += [ ring.radius() ]
            self.__nPhis__           += [ ring.modulesInPhi() ]
            self.__rphistaggers__    += [ ring.phiStagger() ]
            self.__minHalfLengtshX__ += [ ring.module().minHalfX() ]
            self.__maxHalfLengtshX__ += [ ring.module().maxHalfX() ]
            self.__halfLengtshX__    += [ ring.module().halfY() ]
            self.__moduleThickness__ += [ ring.module().thickness( )]
        
    def zPosition(self) :
        return self.__zPosition__

    def rStagger(self) :
        return self.__rStagger__

    def rEnvelope(self) :
        return self.__rEnvelope__

    def rings(self) :
        return self.__rings__
    
    def ringRadii(self) :
        return self.__rRadii__

    def ringPhiModules(self) :
        return self.__nPhis__
        
    def ringModulesPhiStagger(self) :
        return self.__rphistaggers__    
        
    def ringModulesMinHalfX(self) :
        return self.__minHalfLengtshX__
        
    def ringModulesMaxHalfX(self) :
        return self.__maxHalfLengtshX__
        
    def ringModulesHalfY(self) : 
        return self.__halfLengtshX__ 
        
    def ringModulesThickness(self):
        return self.__moduleThickness__  
        
        
""" Barrel type Volume definition """
class BarrelVolume(object):    
    """docstring for BarrelVolume"""
    def __init__(self, barrelLayers):
        super(BarrelVolume, self).__init__()
        self.__barrelLayers__ = barrelLayers
        # process
        self.__layerRadii__            = []
        self.__layerEnvelopeZ__        = []
        self.__layerModulesPhi__       = []
        self.__layerMoudlesTiltPhi__   = []
        self.__layerModulesPositionZ__ = []
        self.__layerModuleStaggerZ__   = []
        self.__layerModuleHalfX__      = []
        self.__layerModuleHalfY__      = []
        self.__layerModuleThickness__  = []
        for layer in self.__barrelLayers__ :
            self.__layerRadii__            += [ layer.radius() ]
            self.__layerEnvelopeZ__        += [ layer.module() ]
            self.__layerModulesPhi__       += [ layer.modulesInPhi() ]
            self.__layerMoudlesTiltPhi__   += [ layer.modulePositionPhi() ]
            self.__layerModulesPositionZ__ += [ layer.modulesInZ() ]
            self.__layerModuleStaggerZ__   += [ layer.modulePositionZ() ]
            self.__layerModuleHalfX__      += [ layer.tiltPhi() ]
            self.__layerModuleHalfY__      += [ layer.overlapZ() ]
            self.__layerModuleThickness__  += [ layer.staggerZ() ]
            
    def layerRadii( self ): 
        return self.__layerRadii__            

    def layerEnvelopesZ( self ): 
        return self.__layerEnvelopeZ__        

    def layerModulesPhi( self ) :
        return self.__layerModulesPhi__       

    def layerModulesTiltPhi( self ) :
        return self.__layerMoudlesTiltPhi__   

    def layerModulesPositionsZ( self ) :
        return self.__layerModulesPositionZ__ 

    def layerModulesStaggerZ( self ) :
        return self.__layerModuleStaggerZ__   

    def layerModulesHalfX( self ) :
        return self.__layerModuleHalfX__      

    def layerModulesHalfY( self ) :
        return self.__layerModuleHalfY__      

    def layerModulesThickness( self ) :
        return self.__layerModuleThickness__  

""" Endcap type Volume definition """
class EndcapVolume(object):
    """docstring for EndcapVolume"""
    def __init__(self, discLayers):
        super(EndcapVolume, self).__init__()
        self.__discLayers__ = discLayers
        # process
        self.__layerPositionsZ__        = []
        self.__layerEnvelopesR__        = []
        self.__layerModulesRadii__      = []
        self.__layerModulesStaggerR__   = []
        self.__layerModulesPhi__        = []
        self.__layerModulesStaggerPhi__ = []
        self.__layerModulesMinHalfX__   = []
        self.__layerModulesMaxHalfX__   = []
        self.__layerModulesHalfY__      = []
        self.__layerModulesThickness__  = []
        for layer in self.__discLayers__ :
            self.__layerPositionsZ__        += [ layer.zPosition() ]
            self.__layerEnvelopesR__        += [ layer.rEnvelope() ]
            self.__layerModulesRadii__      += [ layer.ringRadii() ]
            self.__layerModulesStaggerR__   += [ layer.rStagger()  ]
            self.__layerModulesPhi__        += [ layer.ringPhiModules() ]
            self.__layerModulesStaggerPhi__ += [ layer.ringModulesPhiStagger() ]
            self.__layerModulesMinHalfX__   += [ layer.ringModulesMinHalfX() ]
            self.__layerModulesMaxHalfX__   += [ layer.ringModulesMaxHalfX() ]
            self.__layerModulesHalfY__      += [ layer.ringModulesHalfY() ]
            self.__layerModulesThickness__  += [ layer.ringModulesThickness() ]
        
    def layerPositionsZ(self):                 
        return self.__layerPositionsZ__
    
    def layerEnvelopesR(self):
        return self.__layerEnvelopesR__        
    
    def layerModulesRadii(self): 
        return self.__layerModulesRadii__              
    
    def layerModulesStaggerR(self):
        return self.__layerModulesStaggerR__            
    
    def layerModulesPhi(self):
        return self.__layerModulesPhi__               
    
    def layerModulesStaggerPhi(self):
        return self.__layerModulesStaggerPhi__                     
    
    def layerModulesMinHalfX(self):
        return self.__layerModulesMinHalfX__    
    
    def layerModulesMaxHalfX(self):    
        return self.__layerModulesMaxHalfX__
    
    def layerModulesHalfY(self):
        return self.__layerModulesHalfY__    
           
    def layerModulesThickness(self):
        return self.__layerModulesThickness__
        
           