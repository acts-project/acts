// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ILayerArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TOOLS_ILAYERARRAYCREATOR_H
#define ACTS_TOOLS_ILAYERARRAYCREATOR_H 1

#include <memory>
#include <vector>
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinningType.hpp"

namespace Acts {

class Layer;
/// A std::shared_ptr to a Layer
typedef std::shared_ptr<const Layer> LayerPtr;
/// A BinnedArray to a std::shared_ptr of a layer
typedef BinnedArray<LayerPtr> LayerArray;
/// A vector of std::shared_ptr to layers
typedef std::vector<LayerPtr> LayerVector;

/// @class ILayerArrayCreator
///
/// Interface class ILayerArrayCreators, it inherits from IAlgTool.
///
/// It receives the LayerVector and creaets an array with NaivgationLayer
/// objects
/// filled in between.
class ILayerArrayCreator
{
public:
  /// Virtual destructor
  virtual ~ILayerArrayCreator() = default;

  /// LayerArrayCreator interface method
  ///
  /// @param layers are the layers to be moved into an array
  /// @param min is the minimul value for binning
  /// @param max is the maximum value for binning
  /// @param btype is the binning type
  /// @param bvalue is the value in which the binning should be done
  ///
  /// @return unqiue pointer to a new LayerArray
  virtual std::unique_ptr<const LayerArray>
  layerArray(const LayerVector& layers,
             double             min,
             double             max,
             BinningType        btype  = arbitrary,
             BinningValue       bvalue = binX) const = 0;
};
}  // end of namespace

#endif  // ACTS_TOOLS_ILAYERARRAYCREATOR_H
