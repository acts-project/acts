// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/VolumeBounds.hpp"

#include <cassert>
namespace Acts {
std::ostream& operator<<(std::ostream& sl, const VolumeBounds& vb) {
  return vb.toStream(sl);
}

std::ostream& operator<<(std::ostream& sl, const VolumeBounds::BoundsType& bt) {
  switch (bt) {
    using enum VolumeBounds::BoundsType;
    case eCone:
      sl << "Cone";
      break;
    case eCuboid:
      sl << "Cuboid";
      break;
    case eCutoutCylinder:
      sl << "CutoutCylinder";
      break;
    case eCylinder:
      sl << "Cylinder";
      break;
    case eGenericCuboid:
      sl << "GenericCuboid";
      break;
    case eTrapezoid:
      sl << "Trapezoid";
      break;
    case eDiamond:
      sl << "Diamond";
      break;
    case eOther:
      sl << "Other";
      break;
  }
  return sl;
}

std::vector<OrientedSurface> VolumeBounds::boundarySurfaces(
    const GeometryContext& gctx, Volume& parentVolume) const {
  if (!parentVolume.volumePositioner()) {
    return orientedSurfaces(parentVolume.localToGlobalTransform(gctx));
  }
  std::vector<OrientedSurface> portalSurfaces =
      orientedSurfaces(Acts::Transform3::Identity());

  for (std::size_t faceIdx = 0lu; faceIdx < portalSurfaces.size(); ++faceIdx) {
    std::shared_ptr<RegularSurface>& alignMe = portalSurfaces[faceIdx].surface;
    const Acts::Transform3 internalTrf = alignMe->localToGlobalTransform(gctx);
    alignMe = parentVolume.volumePositioner()->alignWithVolume(
        faceIdx, internalTrf, std::move(alignMe));
    // Ensure that the surface is not destroyed by the client
    if (alignMe == nullptr) {
      throw std::logic_error(
          "boundarySurfaces() - alignWithVolume() must not destroy the "
          "surface");
    }
    // if the surface does not have an associated detector element
    // it cannot move with the volume
    if (alignMe->associatedDetectorElement() == nullptr) {
      throw std::logic_error(
          "boundarySurfaces() - alignWithVolume() is supposed to make the "
          "surface alignable");
    }
    if (alignMe->isSensitive()) {
      throw std::logic_error(
          "boundarySurfaces() - The aligned surface shall not be sensitive");
    }
  }
  return portalSurfaces;
}
}  // namespace Acts
