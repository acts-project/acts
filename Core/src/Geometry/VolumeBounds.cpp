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
    Volume& parentVolume) const {
  if (parentVolume.volumePositioner() == nullptr) {
    return orientedSurfaces(*parentVolume.m_transform);
  }
  std::vector<OrientedSurface> portalSurfaces =
      orientedSurfaces(Acts::Transform3::Identity());
  // It is safe to construct a default geometry context here as the oriented
  // surfaces have their own transform
  const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  for (std::size_t faceIdx = 0lu; faceIdx < portalSurfaces.size(); ++faceIdx) {
    std::shared_ptr<RegularSurface>& surface = portalSurfaces[faceIdx].surface;
    // the localToGlobal from the surface is the transform from the portal
    // frame to the volume center
    const Acts::Transform3 internalTrf = surface->localToGlobalTransform(gctx);
    surface = parentVolume.volumePositioner()->makePortalAlignable(
        faceIdx, internalTrf, std::move(surface));
  }
  return portalSurfaces;
}
}  // namespace Acts
