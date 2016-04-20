///////////////////////////////////////////////////////////////////
// GenericOverlapDescriptor.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_GENERICOVERLAPDESCRIPTOR_H
#define ACTS_DETECTOR_GENERICOVERLAPDESCRIPTOR_H 1

// Core module
#include "ACTS/Utilities/AlgebraDefinitions.h"
// Geometry module
#include "ACTS/Utilities/OverlapDescriptor.h"
#include "ACTS/Utilities/Intersection.h"

namespace Acts {

     class Surface;
    
     /**
     @class GenericOverlapDescriptor
     
     Neighbour & bin member based overlap descriptor

     @author Andreas.Salzburger@cern.ch
    */

    class GenericOverlapDescriptor : public OverlapDescriptor {
      public: 
          
        /**Default constructor - surfaceArray is not given */
        GenericOverlapDescriptor(){}
          
        
        /**Virtual destructor*/
        virtual ~GenericOverlapDescriptor(){}
        
        /**Pseudo-constructor*/
        GenericOverlapDescriptor* clone() const { return new GenericOverlapDescriptor(); }
    
        /** get the compatible surfaces 
            - return : a boolean indicating if an actual intersection had been tried
            - fill vector of intersections
            - primary bin surface : sf
            - position & direction : pos, dir
        */
       virtual  bool reachableSurfaces(std::vector<const Surface*>& cSurfaces, 
                                       const Surface& sf,
                                       const Vector3D& pos,
                                       const Vector3D& dir,
                                       int searchType) const override;
       
    };

}

#endif // ACTS_DETECTOR_GENERICOVERLAPDESCRIPTOR_H
