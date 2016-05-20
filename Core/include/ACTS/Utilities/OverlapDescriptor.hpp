// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// OverlapDescriptor.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_IOVERLAPDESCRIPTOR_H
#define ACTS_DETECTOR_IOVERLAPDESCRIPTOR_H

// Core module
#include "ACTS/Utilities/Intersection.hpp"
#include "ACTS/Utilities/Definitions.hpp"
// STD
#include <vector>

namespace Acts {

     class Surface;

     /**
     @class OverlapDescriptor
     
     BaseClass to be overloaded for describing overlaps 
     and next-by elements for the sub-detector implementations.

     It allows to describe the potentially reachable surfaces

     @author Andreas.Salzburger@cern.ch
    */

    class OverlapDescriptor {
      public: 
        /**Default constructor*/
        OverlapDescriptor(){}
        
        /**Virtual destructor*/
        virtual ~OverlapDescriptor(){}
        
        /**Pseudo-constructor*/
        virtual OverlapDescriptor* clone() const  = 0;
    
        /** get the compatible surfaces 
            - return : a boolean indicating if an actual intersection had been tried
            - fill vector of intersections
            - primary bin surface : sf
            - position & direction : pos, dir
        
           Possible search type given by the Layer :
            1,2 - provide bin surface and registered neighbours and bin mates
            3,4 - provide bin surface and next bin surfaces (if they differ)
            5 - whatever the overlap descriptor returns with this
        
        */
        virtual bool reachableSurfaces(std::vector<const Surface*>& cSurfaces, 
                                       const Surface& sf,
                                       const Vector3D& pos,
                                       const Vector3D& dir,
                                       int searchType) const = 0;
       
    };

}

#endif


