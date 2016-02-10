class DetectorModule(object):
    
    """docstring for DetectorModule: defines a simple detector moduel for generic detetors"""
    def __init__(self, arg):
        super(DetectorModule, minHalfX, maxHalfX, halfY, thickness).__init__()
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


class CylinderLayer(object):
    """docstring for CylinderLayer"""
    def __init__(self, radius, module, numPhi, numZ, tiltPhi, overlapZ, staggerZ ):
        super(CylinderLayer, self).__init__()
        self.__radius__   =   radius 
        self.__module__   =   module 
        self.__numPhi__   =   numPhi 
        self.__numZ__     =   numZ 
        self.__tiltPhi__  =   tiltPhi 
        self.__overlapZ__ =   overlapZ
        self.__staggerZ__ =   staggerZ
        
        
        def radius(self):
            return self.__radius__
        
        def module(self):
            return self.__module__
             
        numPhi 
        numZ 
        tiltPhi 
        overlapZ
        staggerZ