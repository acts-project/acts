// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"

#include <memory>

namespace Acts {

class ConeBounds;
class SurfaceArray;

/// @class ConeLayer
///
/// Class to describe a conical detector layer for tracking, it inherits from
/// both, Layer base class and ConeSurface class
class ConeLayer : virtual public ConeSurface, public Layer {
 public:
  /// Factory for shared layer
  ///
  /// @param transform is the 3D transform that positions the layer in 3D frame
  /// @param cbounds is the conical bound description
  /// @param surfaceArray is the array of sensitive surfaces
  /// @param thickness is the layer thickness along the normal axis
  /// @param ad is the approach descriptor for navigation towards the layer
  /// @param laytyp is the layer type
  ///
  /// @todo change od and ad to unique_ptr
  ///
  /// @return is a shared pointer to a layer
  static std::shared_ptr<ConeLayer> create(
      const Transform3& transform, std::shared_ptr<const ConeBounds> cbounds,
      std::unique_ptr<SurfaceArray> surfaceArray, double thickness = 0.,
      std::unique_ptr<ApproachDescriptor> ad = nullptr,
      LayerType laytyp = active);

  ConeLayer() = delete;
  ConeLayer(const ConeLayer& cla) = delete;
  ~ConeLayer() override = default;
  ConeLayer& operator=(const ConeLayer&) = delete;

  /// Transforms the layer into a Surface representation for extrapolation
  /// @return Const reference to the cone surface representing this layer
  const ConeSurface& surfaceRepresentation() const override;

  /// Non-const version of surface representation access
  /// @return Mutable reference to the cone surface
  ConeSurface& surfaceRepresentation() override;

 protected:
  /// Private constructor with arguments
  ///
  /// @param transform is the 3D transform that positions the layer in 3D frame
  /// @param cbounds is the conical bound description
  /// @param surfaceArray is the array of sensitive surfaces
  /// @param thickness is the layer thickness along the normal axis
  /// @param ade is the approach descriptor for navigation towards the layer
  /// @param laytyp is the layer type
  ///
  /// @todo change od and ad to unique_ptr
  ConeLayer(const Transform3& transform,
            std::shared_ptr<const ConeBounds> cbounds,
            std::unique_ptr<SurfaceArray> surfaceArray, double thickness = 0.,
            std::unique_ptr<ApproachDescriptor> ade = nullptr,
            LayerType laytyp = active);

  /// Private copy constructor with shift, called by create(args*)
  ///
  /// @param cla is the source clone layer
  /// @param shift is the additional shift applied after copying
  ConeLayer(const ConeLayer& cla, const Transform3& shift);
};

}  // namespace Acts
