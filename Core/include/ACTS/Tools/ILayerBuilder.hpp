// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ILayerBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ILAYERBUILDER_H
#define ACTS_GEOMETRYINTERFACES_ILAYERBUILDER_H 1

// STL
#include <vector>
#include <string>
#include <memory>

namespace Acts {

  class Layer;
  typedef std::shared_ptr<const Layer> LayerPtr;
  typedef std::vector< LayerPtr > LayerVector;

  /** @class ILayerBuilder
  
    Interface class for ILayerBuilders in a typical 
    | EC- | Central | EC+ | 
    detector setup.
      
    */
  class ILayerBuilder
  { 
    public:
      /**Virtual destructor*/
      virtual ~ILayerBuilder(){}

      /** LayerBuilder interface method - returning the layers at negative side */
      virtual const LayerVector negativeLayers() const = 0;
      
      /** LayerBuilder interface method - returning the central layers */
      virtual const LayerVector centralLayers() const = 0;
      
      /** LayerBuilder interface method - returning the layers at negative side */
      virtual const LayerVector positiveLayers() const = 0; 

      /** Name identification */
      virtual const std::string& identification() const = 0;
             
  };

} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_ILAYERBUILDER_H
