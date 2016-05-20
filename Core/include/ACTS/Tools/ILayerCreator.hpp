// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ILayerCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ILAYERCREATOR_H
#define ACTS_GEOMETRYINTERFACES_ILAYERCREATOR_H 1

// STL
#include <vector>
#include <string>
#include <memory>

namespace Acts {

  class Surface;
  class Layer;
  typedef std::shared_ptr<const Layer> LayerPtr;

  /** @class ILayerCreator
  
    Interface class for LayerCreator from DetectorElements
      
    */
  class ILayerCreator {
    
    public:
      /**Virtual destructor*/
      virtual ~ILayerCreator(){}

      /** ILayerCreator interface method - returning a cylindrical layer */
      virtual LayerPtr cylinderLayer(const std::vector<const Surface*>& surfaces,
                                     double envelopeR, double evelopeZ,
                                     size_t binsPhi, size_t binsZ) const = 0; 
    
      /** ILayerCreator interface method - returning a disc layer */
      virtual LayerPtr discLayer(const std::vector<const Surface*>& surfaces,
                                 double envelopeMinR, double envelopeMaxR, double envelopeZ,
                                 size_t binsR, size_t binsPhi,
                                 const std::vector<double>& rBoundaries = {}) const = 0; 
    
      /** ILayerCreator interface method - returning a plane layer */
      virtual LayerPtr planeLayer(const std::vector<const Surface*>& surfaces,
                                  double envelopeXY, double envelopeZ,
                                  size_t binsX, size_t binsY) const = 0; 
             
  };

} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_ILAYERCREATOR_H
