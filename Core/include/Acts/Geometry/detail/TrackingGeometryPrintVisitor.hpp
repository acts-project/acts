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

#include <sstream>

namespace Acts::detail {
/// @brief  Visitor class which prints the content of the tracking geometry (volumes, surfaces, portals) in an
///         indented list. The indentation reflects the position of a volume
///         inside the overall tracking geometry
class TrackingGeometryPrintVisitor : public Acts::TrackingGeometryVisitor {
 public:
  /// @brief Constructor
  /// @param gctx: Reference to the geometry context needed to align the surfaces inside the detector
  /// @param indentation: How many spaces indicate a new step inside the volume hierarchy.
  explicit TrackingGeometryPrintVisitor(const Acts::GeometryContext& gctx,
                                        std::size_t indentation = 4);
  ~TrackingGeometryPrintVisitor() override;

  /// @brief Visit a tracking volume in the geometry
  /// @param volume The tracking volume being visited
  /// @note Called for each volume in the geometry hierarchy during traversal
  void visitVolume(const TrackingVolume& volume) override;

  /// @brief Visit and potentially modify a portal
  /// @param portal The portal being visited
  /// @note Called for each portal encountered during geometry traversal
  void visitPortal(const Portal& portal) override;

  /// @brief Visit and potentially modify a surface
  /// @param surface The surface being visited
  /// @note Called for each surface encountered during geometry traversal
  void visitSurface(const Surface& surface) override;
  /// @brief Return the reference to the underlying stream containing the printout
  std::stringstream& stream();
  /// @brief Return the reference to the underlying stream containing the printout
  const std::stringstream& stream() const;

 private:
  /// @brief Calculates the depth of the volume in the hierarchy tree
  void updateDepth(const TrackingVolume& trkVol);
  /// @brief Returns the position of the volume in the parent tree
  std::size_t volNumber(const TrackingVolume& trkVol) const;
  /// @brief Geometry context to fetch the center positions of the surface
  Acts::GeometryContext m_gctx;
  /// @brief Stream into which all print-out is piped
  std::stringstream m_printStream;
  /// @brief Current depth in the volume hierarchy
  std::size_t m_currentDepth{0};
  /// @brief Indentiation size
  std::size_t m_indentation{4};
};
}  // namespace Acts::detail
