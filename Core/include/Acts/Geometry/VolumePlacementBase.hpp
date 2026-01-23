// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/detail/PortalPlacement.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"

#include <memory>
#include <vector>

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
  /// @brief Default constructor
  VolumePlacementBase();
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

  virtual bool portalTransformCached(const std::size_t portalIdx) const = 0;

  virtual void cachePortalTransform(const std::size_t portalIdx,
                                    const Transform3& portalLocToGlob) = 0;

  virtual const Transform3& portalLocalToGlobal(
      const GeometryContext& gctx, const std::size_t portalIdx) const = 0;

  std::shared_ptr<RegularSurface> makePortalAlignable(
      const std::size_t portalIdx, const Transform3& portalToVolTrf,
      std::shared_ptr<RegularSurface>&& surface);

 private:
  std::vector<std::unique_ptr<detail::PortalPlacement>> m_portalPlacements{};
};
}  // namespace Acts
