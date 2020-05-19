// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ApproachDescriptor.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

class Surface;
class Layer;
class BoundaryCheck;

/// @class ApproachDescriptor
///
/// Virtual base class to decide and return which approaching surface to be
/// taken, the surfaces are std::shared_ptr, as they can be the boundary
/// surfaces of the representingVolume of the Layer
class ApproachDescriptor {
 public:
  ApproachDescriptor() = default;

  /// Virtual destructor
  virtual ~ApproachDescriptor() = default;

  /// @brief Register Layer
  /// Links the layer to the approach surfaces
  ///
  /// @param lay is the layer to be assigned
  virtual void registerLayer(const Layer& lay) = 0;

  /// @brief Get the surface on approach
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the position from start of the search
  /// @param direction is the direction at the start of the search
  /// @param bcheck is the boundary check directive
  ///
  /// @return is a surface intersection
  virtual ObjectIntersection<Surface> approachSurface(
      const GeometryContext& gctx, const Vector3D& position,
      const Vector3D& direction, const BoundaryCheck& bcheck) const = 0;

  /// Get all the contained surfaces
  /// @return all contained surfaces of this approach descriptor
  virtual const std::vector<const Surface*>& containedSurfaces() const = 0;

  /// Non-const version
  virtual std::vector<const Surface*>& containedSurfaces() = 0;
};

}  // namespace Acts
