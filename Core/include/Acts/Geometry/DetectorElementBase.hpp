// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/SurfacePlacementBase.hpp"

namespace Acts {

/// This class is the base of all detector elements that are usable by ACTS.
/// All experiment-specific detector element classes are expected to inherit
/// from it.
///
/// @deprecated: This class is deprecated in favour of SurfacePlacementBase
/// @remark It is possible toe replace this base class by defining a
///         `ACTS_DETECTOR_ELEMENT_BASE_REPLACEMENT` pre-processor replacement.
///         If found, @ref SurfacePlacementBase.hpp will instead include that file.
[[deprecated(
    "This class is deprecated in favour of SurfacePlacementBase")]] class
    DetectorElementBase : public SurfacePlacementBase {
 public:
  DetectorElementBase() = default;
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
      const GeometryContext& gctx) const override {
    ACTS_PUSH_IGNORE_DEPRECATED()
    return transform(gctx);
    ACTS_POP_IGNORE_DEPRECATED()
  }
  /// Returns the thickness of the module
  /// @return double that indicates the thickness of the module
  virtual double thickness() const = 0;
  /// Returns whether the detector element corresponds to a sensitive
  /// surface on which measurements are expressed
  virtual bool isSensitive() const override { return true; }
};

}  // namespace Acts

#endif
