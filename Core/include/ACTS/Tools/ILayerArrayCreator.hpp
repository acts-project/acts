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
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/BinningType.hpp"

namespace Acts {

class Layer;
typedef std::shared_ptr<const Layer> LayerPtr;

/// @typedef LayerArray and LayerVector
typedef BinnedArray<LayerPtr> LayerArray;
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
  /// @param layers are the layers to be moved into an array
  /// @min is the minimul value for binning
  /// @max is the maximum value for binning
  /// @btype is the binning type
  /// @bvalue is the value in which the binning should be done
  virtual std::unique_ptr<const LayerArray>
  layerArray(const LayerVector& layers,
             double             min,
             double             max,
             BinningType        btype  = arbitrary,
             BinningValue       bvalue = binX) const = 0;
};
}  // end of namespace

#endif  // ACTS_TOOLS_ILAYERARRAYCREATOR_H
