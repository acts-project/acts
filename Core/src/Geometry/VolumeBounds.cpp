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
  if (!parentVolume.isAlignable()) {
    return orientedSurfaces(*parentVolume.m_transform);
  }
  return parentVolume.volumePlacement()->makePortalsAlignable(
      orientedSurfaces(Acts::Transform3::Identity()));
}
}  // namespace Acts
