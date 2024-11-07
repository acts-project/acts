// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinnedArray.hpp"

#include <memory>
#include <vector>

namespace Acts {

class Layer;
/// A std::shared_ptr to a Layer
using LayerPtr = std::shared_ptr<const Layer>;
/// A BinnedArray to a std::shared_ptr of a layer
using LayerArray = BinnedArray<LayerPtr>;
/// A vector of std::shared_ptr to layers
using LayerVector = std::vector<LayerPtr>;

/// @class ILayerArrayCreator
///
/// Interface class ILayerArrayCreators.
///
/// It receives the LayerVector and creaets an array with NaivgationLayer
/// objects filled in between.
class ILayerArrayCreator {
 public:
  /// Virtual destructor
  virtual ~ILayerArrayCreator() = default;

  /// LayerArrayCreator interface method
  ///
  /// @param gctx is the geometry context for witch the volume is built
  /// @param layers are the layers to be moved into an array
  /// @param min is the minimul value for binning
  /// @param max is the maximum value for binning
  /// @param bType is the type of the ordering: Equidistant, Variable
  /// @param bValue is the value in which the layer ordering is done
  ///
  /// @return unique pointer to a new LayerArray
  virtual std::unique_ptr<const LayerArray> layerArray(
      const GeometryContext& gctx, const LayerVector& layers, double min,
      double max, AxisType bType = AxisType::Variable,
      AxisDirection bValue = AxisDirection::AxisX) const = 0;
};
}  // namespace Acts
