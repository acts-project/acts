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
#include "Acts/Utilities/BinUtility.hpp"
#include "ActsFatras/Digitization/PlanarSurfaceDrift.hpp"
#include "ActsFatras/Digitization/PlanarSurfaceMask.hpp"
#include "ActsFatras/Digitization/Segmentizer.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <numeric>

namespace Acts {
class Surface;
}

namespace ActsFatras {

/// @brief Class that ties the digitization modules together and produces the channels
class Channelizer {
  PlanarSurfaceDrift m_surfaceDrift;
  PlanarSurfaceMask m_surfaceMask;
  Segmentizer m_segmentizer;

 public:
  /// Do the geometric channelizing
  ///
  /// @param hit The hit we want to channelize
  /// @param surface the surface on which the hit is
  /// @param gctx the Geometry context
  /// @param driftDir the drift direction
  /// @param segmentation the segmentation of the surface
  /// @param thickness the thickness of the surface
  /// @param minRelPerpDrift minimum relative perpendicular drift (to avoid numerical instability)
  ///
  /// @return the list of channels
  Acts::Result<std::vector<Segmentizer::ChannelSegment>> channelize(
      const Hit& hit, const Acts::Surface& surface,
      const Acts::GeometryContext& gctx, const Acts::Vector3& driftDir,
      const Acts::BinUtility& segmentation, double thickness,
      double minRelPerpDrift = 0.001) const;
};

}  // namespace ActsFatras
