// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

namespace Acts {
class Surface;
/// @brief The `SurfacePlacementBase` is an API proxy to model the dynamic
///        movement of the `Acts::Surface` representing e.g. the experiment's
///        sensor planes.
///
///        If the user wants to let the sensor surface float, he has to
///        implement an experiment-specifc class inheriting from
///        SurfacePlacementBase and to instantiate the surface with such an
///        instance.
///
///        In every call for the Surface's current position the GeometryContext
///        which is also contextualized by the experiment is piped through the
///        entire call stack until its piped back to the experiment specific
///        implementation via the call of the
///        `SurfacePlacementBase::localToGlobalTransform`. There the user needs
///        to unpack the GeometryContext, look-up the appropriate cached
///        transform and return it back to the Acts library
class SurfacePlacementBase {
 public:
  /// @brief Virtual default constructor
  virtual ~SurfacePlacementBase() = default;

  /// @brief Return the transform for the Element proxy mechanism
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return reference to the transform to switch from the element's
  ///         coordinates to the experiment's global coordinate system
  virtual const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const = 0;

  /// @brief Get a reference to the surface that is associated
  ///        with this placement element.
  /// @note It is expected that the surface returned will have it's @ref
  ///       Acts::Surface::surfacePlacement method return a pointer to
  ///       this object.
  /// @return Reference to a surface that represents this detector element
  virtual const Surface& surface() const = 0;

  /// @copydoc surface
  /// @return Reference to a surface that represents this detector element
  virtual Surface& surface() = 0;

  /// @brief Returns whether the placement corresponds to a surface on which
  ///        the measurements from the experiment are represented, i.e. it is
  //         a detector surface
  /// @return True if this is a sensitive surface
  virtual bool isSensitive() const = 0;
};
}  // namespace Acts
