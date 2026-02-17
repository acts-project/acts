// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometryVisitor.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts::detail {
/// TrackingGeometry visitor that connects the portal surfaces of the
/// alignable volume with the volume's alignment system (cf.
/// VolumePlacementBase) The portal surfaces of the alignable volumes are
/// collected based on the geometry identifier and then parsed to the
/// VolumePlacementBase which takes over and creates the backend caches
/// for their aligned transform
class AlignablePortalVisitor final : public TrackingGeometryMutableVisitor {
 public:
  /// Standard constructor taking the geometry context with the current
  /// alignment
  /// @param gctx : Geometry context containing the current alignment configuration
  /// @param logger: Reference to the logger of the object calling the visitor
  AlignablePortalVisitor(const GeometryContext& gctx, const Logger& logger);
  ///  Enables the visitor's inspection mode. The transforms of the portals
  ///  before and after the alignment procedure are required to be
  ///  equivalent Warnings are printed if the tracking volume contains
  ///  un-aligneable subvolumes
  void enableInspection();
  /// Visit and potentially modify a tracking volume
  /// @param volume The tracking volume being visited
  /// @note Called for each volume in the geometry hierarchy during traversal
  void visitVolume(TrackingVolume& volume) override;

 private:
  ///  The construction geometry context connected to the experiment's alignment
  const GeometryContext m_gctx{GeometryContext::dangerouslyDefaultConstruct()};
  ///  Pointer to the held logger object
  const Logger& m_logger;
  ///  Inspect the portal transforms
  bool m_doInspect{false};
  ///  Return the reference to the logger object
  const Logger& logger() const;
};

}  // namespace Acts::detail
