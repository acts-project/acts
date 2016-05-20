// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ApproachDescriptor.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_APPROACHDESCRIPTOR_H
#define ACTS_GEOMETRYUTILS_APPROACHDESCRIPTOR_H 1

// Core module
#include "ACTS/Utilities/Intersection.hpp"
#include "Definitions.hpp"

namespace Acts {

    class Surface;
    class Layer;
    class BoundaryCheck;
    class ICompatibilityEstimator;
    
    typedef ObjectIntersection<Surface> SurfaceIntersection;
         
     /**
     @class ApproachDescriptor
     
     Virtual base class to decide and return which approaching surface to be taken,
     the surfaces are std::shared_ptr, as they can be the boundary surfaces of the 
     representingVolume of the Layer
     
     @author Andreas.Salzburger@cern.ch 
     */

    class ApproachDescriptor {
      public: 
        /** Default constructor */
        ApproachDescriptor(){}
                
        /** Virtual destructor */
        virtual ~ApproachDescriptor(){}

        /** register Layer */
        virtual void registerLayer(const Layer& lay) = 0;
        
        /** get the surface on approach - nullptr means that there's no surface on approach */
        virtual const SurfaceIntersection approachSurface(const Vector3D& pos, 
                                                          const Vector3D& dir, 
                                                          const BoundaryCheck& bchk,
                                                          const ICompatibilityEstimator* ice = nullptr) const = 0;

        /** return all contained surfaces of this approach descriptor */
        virtual const std::vector< const Surface* >& containedSurfaces() const = 0;
        
    };
}

#endif // ACTS_GEOMETRYUTILS_APPROACHDESCRIPTOR_H
