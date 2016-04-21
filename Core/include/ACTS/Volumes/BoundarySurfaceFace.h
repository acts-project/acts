///////////////////////////////////////////////////////////////////
// BoundarySurfaceFace.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_BOUNDARYSURFACEFACE_H
#define ACTS_VOLUMES_BOUNDARYSURFACEFACE_H 1
  
namespace Acts {  
  
  /**
    @enum BoundarySurfaceFace
    
    Enum to describe the position of the BoundarySurface 
    respectively to the frame orientatin of the volume,
    this is mainly ment for code readability.
    
    The different numeration sequences can be found in the
    documentation of the actual VolumeBounds implementations.
    
    The order of faces is chosen to follow - as much as 
    possible - a cycular structure.
    
    @author Andreas.Salzburger@cern.ch 
    
    */

    enum BoundarySurfaceFace {  negativeFaceXY         = 0,  
                                positiveFaceXY         = 1,
                                negativeFaceYZ         = 2,
                                positiveFaceYZ         = 3,
                                negativeFaceZX         = 4,
                                positiveFaceZX         = 5,
                                cylinderCover          = 2,
                                tubeInnerCover         = 3,
                                tubeOuterCover         = 2,
                                tubeSectorNegativePhi  = 4,
                                tubeSectorPositivePhi  = 5,
                                tubeSectorInnerCover   = 3,
                                tubeSectorOuterCover   = 2,
                                trapezoidFaceAlpha     = 2,
                                trapezoidFaceBeta      = 3,
                                index0                 = 0,     
                                index1                 = 1,     
                                index2                 = 2,     
                                index3                 = 3,     
                                index4                 = 4,     
                                index5                 = 5,     
                                index6                 = 6,     
                                index7                 = 7,     
                                index8                 = 8,     
                                index9                 = 9,     
                                index10                = 10,     
                                index11                = 11,
                                undefinedFace          = 99

    };
     
    /** @enum to specify the inside/outside
        this is with respect to the normal vector
     */
    enum BoundaryOrientation {
        insideVolume  = -1,
        outsideVolume = 1
    };
     
}

#endif // ACTS_VOLUMES_BOUNDARYSURFACEFACE_H

