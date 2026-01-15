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

#ifdef ACTS_DETECTOR_ELEMENT_BASE_REPLACEMENT
#include ACTS_DETECTOR_ELEMENT_BASE_REPLACEMENT
#else

namespace Acts {

class Surface;

/// This class is the base of all detector elements that are usable by ACTS.
/// All experiment-specific detector element classes are expected to inherit
/// from it.
///
/// @remark It is possible toe replace this base class by defining a
///         `ACTS_DETECTOR_ELEMENT_BASE_REPLACEMENT` pre-processor replacement.
///         If found, @ref DetectorElementBase.hpp will instead include that file.
class DetectorElementBase {
 public:
  DetectorElementBase() = default;
  virtual ~DetectorElementBase() = default;

  /// Return the transform for the Element proxy mechanism
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return reference to the transform of this detector element
  [[deprecated(
      "Please use localToGlobalTransform(const GeometryContext& gctx) "
      "instead")]]
  virtual const Transform3& transform(const GeometryContext& gctx) const {
    return localToGlobalTransform(gctx);
  }

  /// Return the transform for the Element proxy mechanism
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return reference to the transform to switch from the element's
  ///         coordinates to the experiment's global coordinate system
  virtual const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const {
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    return transform(gctx);
#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
  }

  /// Get a reference to the surface that is associated with this detector
  /// element.
  /// @note It is expected that the surface returned will have it's @ref
  ///       Acts::Surface::associatedDetectorElement method return a pointer to
  ///       this object.
  /// @return Reference to a surface that represents this detector element
  virtual const Surface& surface() const = 0;

  /// @copydoc surface
  /// @return Reference to a surface that represents this detector element
  virtual Surface& surface() = 0;

  /// Returns the thickness of the module
  /// @return double that indicates the thickness of the module
  virtual double thickness() const = 0;
};

}  // namespace Acts

#endif
