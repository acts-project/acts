// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_CYLINDERLAYER_H
#define ACTS_DETECTOR_CYLINDERLAYER_H

class MsgStream;

#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Definitions.hpp"
// STL sorting
#include <algorithm>

namespace Acts {

class CylinderBounds;
class SurfaceMaterial;
class OverlapDescriptor;
class ApproachDescriptor;

  /**
   @class CylinderLayer
   
   Class to describe a cylindrical detector layer for tracking, it inhertis from both, 
   Layer base class and CylinderSurface class
       
  */

  class CylinderLayer : public CylinderSurface, public Layer {
                   
    public:
      /**create a shared, fully deployed CylinderLayer */
      static LayerPtr create(std::shared_ptr<Transform3D> transform,
                             std::shared_ptr<const CylinderBounds> cbounds,
                             std::unique_ptr<SurfaceArray> surfaceArray = nullptr,
                             double thickness = 0.,
                             OverlapDescriptor* od  = nullptr,
                             ApproachDescriptor* ad = nullptr,
                             int laytyp=int(Acts::passive))
      { return LayerPtr(new CylinderLayer(transform, cbounds, std::move(surfaceArray), thickness, od, ad, laytyp));}
        
      /** Copy constructor with shift - will not copy anything of the static next/previous environment*/
      static LayerPtr create(const CylinderLayer& cla, const Transform3D& shift) 
      { return LayerPtr(new CylinderLayer(cla, shift)); }                                                           
      
      /** Clone with a shift - only cloning that is allowed */
      LayerPtr cloneWithShift(const Transform3D& shift) const override
      { return CylinderLayer::create(*this,shift); }
      
      /** Copy constructor - forbidden, create a new one if you need */
      CylinderLayer(const CylinderLayer& cla) = delete;
       
      /** Assignment operator for CylinderLayers - forbidden, create a new one */
      CylinderLayer& operator=(const CylinderLayer&) = delete;
                    
      /** Destructor*/
      virtual ~CylinderLayer(){}
              
      /** Transforms the layer into a Surface representation for global positioning & navigation */
      const CylinderSurface& surfaceRepresentation() const override;
                                        
     private:   
       /** build approach surfaces */
       void buildApproachDescriptor() const;
       
     protected:   
          
       /** Default Constructor*/
       CylinderLayer(){}
                      
       /** Fully deployed CylinderLayer constructor */
       CylinderLayer(std::shared_ptr<Transform3D> transform,
                     std::shared_ptr<const CylinderBounds> cbounds,
                     std::unique_ptr<SurfaceArray> surfaceArray = nullptr,
                     double thickness = 0.,
                     OverlapDescriptor* od  = nullptr,
                     ApproachDescriptor* ad = nullptr,
                     int laytyp=int(Acts::passive));
       
       /** Copy constructor with shift - will not copy anything of the static next/previous environment*/
       CylinderLayer(const CylinderLayer& cla, const Transform3D& tr);
       
  };

} // end of namespace

#endif // TRKGEOMETY_CYLINDERLAYER_H
