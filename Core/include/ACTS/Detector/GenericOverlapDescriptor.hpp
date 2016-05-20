// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GenericOverlapDescriptor.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_GENERICOVERLAPDESCRIPTOR_H
#define ACTS_DETECTOR_GENERICOVERLAPDESCRIPTOR_H 1

// Core module
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/OverlapDescriptor.hpp"
#include "ACTS/Utilities/Intersection.hpp"

namespace Acts {

     class Surface;
    
     /**
     @class GenericOverlapDescriptor
     
     Neighbour & bin member based overlap descriptor

    */

    class GenericOverlapDescriptor : public OverlapDescriptor {
      public: 
          
        /**Default constructor - surfaceArray is not given */
        GenericOverlapDescriptor(){}
          
        
        /**Virtual destructor*/
        virtual ~GenericOverlapDescriptor(){}
        
        /**Pseudo-constructor*/
        GenericOverlapDescriptor* clone() const override { return new GenericOverlapDescriptor(); }
    
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
