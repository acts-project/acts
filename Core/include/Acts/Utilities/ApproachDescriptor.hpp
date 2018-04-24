// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ApproachDescriptor.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_UTILITIES_APPROACHDESCRIPTOR_H
#define ACTS_UTILITIES_APPROACHDESCRIPTOR_H 1

#include <vector>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

class Surface;
class Layer;
class BoundaryCheck;
class ICompatibilityEstimator;

typedef ObjectIntersection<Surface> SurfaceIntersection;

/// @class ApproachDescriptor
///
/// Virtual base class to decide and return which approaching surface to be
/// taken,
/// the surfaces are std::shared_ptr, as they can be the boundary surfaces of
/// the
/// representingVolume of the Layer
///
///
class ApproachDescriptor
{
public:
  /// Default constructor
  ApproachDescriptor() {}
  /// Virtual destructor
  virtual ~ApproachDescriptor() {}
  /// Register Layer
  /// this gives the approach surfaces the link to the layer
  ///
  /// @param lay is the layer to be assigned
  virtual void
  registerLayer(const Layer& lay)
      = 0;

  // Get the surface on approach
  // - nullptr means that there's no surface on approach
  ///
  /// @param pos is the position from start of the search
  /// @param dir is the direction
  /// @param bcheck is the boundary check directive
  /// @param ice is a (future) compatibility estimator
  ///
  /// @return is a surface isntersection
  virtual const SurfaceIntersection
  approachSurface(const Vector3D&                pos,
                  const Vector3D&                dir,
                  const BoundaryCheck&           bcheck,
                  const ICompatibilityEstimator* ice = nullptr) const = 0;

  /// Tet to all the contained surfaces
  /// @return all contained surfaces of this approach descriptor
  virtual const std::vector<const Surface*>&
  containedSurfaces() const = 0;

  /// Non-const version
  virtual std::vector<const Surface*>&
  containedSurfaces()
      = 0;
};
}

#endif  // ACTS_GEOMETRYUTILS_APPROACHDESCRIPTOR_H
