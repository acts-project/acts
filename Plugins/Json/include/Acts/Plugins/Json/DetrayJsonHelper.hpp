// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <string>
#include <tuple>

namespace Acts {

namespace DetrayJsonHelper {

/// @brief Helper function to switch keys from ACTS to detray
///
/// DETRAY types @todo change to detray imports when available
///    annulus2 = 0u,
///    cuboid3 = 1u,
///    cylinder2 = 2u,
///    cylinder3 = 3u,
///    portal_cylinder2 = 4u,
///    rectangle2 = 5u,
///    ring2 = 6u,
///    trapezoid2 = 7u,
///    cell_wire = 8u,
///    straw_wire = 9u,
///    single1 = 10u,
///    single2 = 11u,
///    single3 = 12u,
///    unknown = 13u
///
/// @param sBounds is the surface bounds type
/// @param portal is the flag for conversion into detray portal format
///
/// @return type and value array in detray format
inline static std::tuple<unsigned int, std::vector<ActsScalar>> maskFromBounds(
    const Acts::SurfaceBounds& sBounds, bool portal = false) {
  auto bType = sBounds.type();
  auto bValues = sBounds.values();
  // Return value
  unsigned int type = 13u;
  std::vector<double> boundaries = bValues;
  // Special treatment for some portals
  if (portal and bType == SurfaceBounds::BoundsType::eCylinder) {
    boundaries = {bValues[0u], -bValues[1u], bValues[1u]};
    type = 4u;
  } else {
    switch (bType) {
      case SurfaceBounds::BoundsType::eAnnulus: {
        type = 0u;
      } break;
      case SurfaceBounds::BoundsType::eRectangle: {
        type = 5u;
        // ACTS: eMinX = 0, eMinY = 1, eMaxX = 2, eMaxY = 3,
        // detray: e_half_x, e_half_y
        boundaries = {0.5 * (bValues[2] - bValues[0]),
                      0.5 * (bValues[3] - bValues[1])};
      } break;
      case SurfaceBounds::BoundsType::eCylinder: {
        boundaries = {bValues[0u], -bValues[1u], bValues[1u]};
        type = 2u;
      } break;
      case SurfaceBounds::BoundsType::eTrapezoid: {
        type = 7u;
        boundaries = {bValues[0u], bValues[1u], bValues[2u],
                      1 / (2 * bValues[2u])};
      } break;
      case SurfaceBounds::BoundsType::eDisc: {
        boundaries = {bValues[0u], bValues[1u]};
        type = 6u;
      } break;
      default:
        break;
    }
  }
  return std::tie(type, boundaries);
}

/// @brief add volume link
///
/// @param jSurface [in,out] is the json object to be patched
/// @param vLink is the volume link to be added
inline static void addVolumeLink(nlohmann::json& jSurface, int vLink) {
  jSurface["volume_link"] = vLink;
}

}  // namespace DetrayJsonHelper
}  // namespace Acts
