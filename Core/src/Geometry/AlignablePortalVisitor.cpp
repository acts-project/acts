// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/detail/AlignablePortalVisitor.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Geometry/detail/TrackingGeometryPrintVisitor.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <cassert>

namespace {
/// @brief checks whether the two transforms are the same
/// @param a: First transform to check
/// @param b: Second transform to check
inline bool isSame(const Acts::Transform3& a, const Acts::Transform3& b) {
  const Acts::Transform3 c = a * b.inverse();
  if (c.translation().norm() > Acts::s_onSurfaceTolerance) {
    return false;
  }
  for (std::size_t d = 0; d < 3; ++d) {
    const Acts::Vector3 e = Acts::Vector3::Unit(d);
    if (std::abs(e.dot(c * e) - 1.) > Acts::s_onSurfaceTolerance) {
      return false;
    }
  }
  return true;
}
}  // namespace

namespace Acts::detail {

AlignablePortalVisitor::AlignablePortalVisitor(const GeometryContext& gctx,
                                               const Logger& logger)
    : m_gctx{gctx}, m_logger{logger} {}

void AlignablePortalVisitor::visitVolume(TrackingVolume& volume) {
  if (!volume.isAlignable()) {
    ACTS_VERBOSE("AlignablePortalVisitor() - the volume "
                 << volume.geometryId() << ", " << volume.volumeBounds()
                 << " is not alignable");
    return;
  }
  ACTS_DEBUG("AlignablePortalVisitor() - Found alignable volume "
             << volume.geometryId() << ", " << volume.volumeBounds());
  std::vector<std::shared_ptr<RegularSurface>> alignable{};
  for (Portal& portal : volume.portals()) {
    ACTS_DEBUG("AlignablePortalVisitor() - Check wheher portal "
               << portal.surface().geometryId() << " is a boundary surface");
    if (portal.surface().geometryId().withBoundary(0) != volume.geometryId()) {
      continue;
    }
    alignable.push_back(
        dynamic_pointer_cast<RegularSurface>(portal.surface().getSharedPtr()));
    assert(alignable.back() != nullptr);
  }
  ACTS_DEBUG("AlignablePortalVisitor() - Associate "
             << alignable.size() << " porttals with the volume.");
  std::vector<Transform3> portalTrfBefore{};
  // If the visitor is in inspection mode, cache the transforms of the
  // unaligned portals
  if (m_doInspect) {
    std::ranges::transform(alignable, std::back_inserter(portalTrfBefore),
                           [&](const std::shared_ptr<RegularSurface>& surface) {
                             return surface->localToGlobalTransform(m_gctx);
                           });
  }
  volume.volumePlacement()->makePortalsAlignable(m_gctx, alignable);

  if (!m_doInspect) {
    return;
  }
  // Ensure that the portal remain where they were supposed to be
  for (std::size_t p = 0; p < alignable.size(); ++p) {
    if (!isSame(alignable[p]->localToGlobalTransform(m_gctx),
                portalTrfBefore[p])) {
      ACTS_ERROR("AlignablePortalVisitor() - the "
                 << p << "-the portal transforms\n --- aligned:   "
                 << toString(alignable[p]->localToGlobalTransform(m_gctx))
                 << "\n --- unaligned: " << toString(portalTrfBefore[p])
                 << " are not the same");
      throw std::runtime_error(
          "AlignablePortalVisitor() - The portal alignment failed");
    }
  }
  // Next check whether the children volumes are also alignable. Don't throw an
  // exception as it's the user's choice to declare the volumes alignable
  if (std::ranges::any_of(volume.volumes(), [](const Volume& childVol) {
        return !childVol.isAlignable();
      })) {
    std::stringstream sstr{};
    TrackingGeometryPrintVisitor printVisitor{m_gctx};
    volume.apply(printVisitor);
    ACTS_WARNING(
        "AlignablePortalVisitor() - The volume has non alignable children: \n"
        << printVisitor.stream().str());
  }
}

const Logger& AlignablePortalVisitor::logger() const {
  return m_logger;
}

void AlignablePortalVisitor::enableInspection() {
  m_doInspect = true;
}

}  // namespace Acts::detail
