// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"

namespace Acts {
class GeometryContext;

/// @brief Interface class to define the backend for alignable volumes
///        that move coherently with the snesitive surfaces inside
///        The interface provides the transform from local -> global
///        coordinates as well as its way back.
///        To move the oriented surfaces along with the volume itself
///        the interface needs to provide a factory that creates the
///        detector elements that are then passed to the oriented surfaces
class VolumePlacementBase {
 public:
  /// @brief Default destructor
  virtual ~VolumePlacementBase() = default;
  /// @brief Returns the transformation from the local volume coordinates to
  ///        the experiment's global coordinate system
  /// @param gctx The current geometry context object, e.g. alignment
  virtual const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const = 0;
  /// @brief Returns the transformation from the experiment's global frame to the
  ///        local volume coordinate system
  /// @param gctx The current geometry context object, e.g. alignment
  virtual const Transform3& globalToLocalTransform(
      const GeometryContext& gctx) const = 0;

  /// @brief Interface method that connects the portal surface to the volume's
  ///        alignment. For the given volume surface face and surface, the
  ///        client code is expected to create an DetectorElementBase and to
  ///        parse it to the surface via
  ///             portalSurface->assignDetectorElement(detElement)
  ///        and then to return back the portal surface again.
  ///        The ownership is transferred to the client as multiple calls on the
  ///        same face are not excluded.
  /// @param faceIdx: Index of the oriented surface
  /// @param internalTrf: Transform aligning the portal within the volume
  virtual std::shared_ptr<RegularSurface> alignWithVolume(
      const std::size_t faceIdx, const Transform3& internalTrf,
      std::shared_ptr<RegularSurface>&& portalSurface) = 0;
};
}  // namespace Acts
