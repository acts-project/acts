// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <memory>
#include <vector>

namespace Acts {

/// @brief  This is a helper factory to create surface objects and
/// and/or bounds from principle parameters
///
/// This can be used e.g. when reading perisitified data and create
/// surfaces from that.
namespace SurfaceFactory {

/// @brief Create a surface object from a transform, a type
/// and bounds parameters and bounds type
///
/// @param transform the transform positions the surface in 3D space
/// @param sType the surface type indicates which surface should be constructed
/// @param bType the bounds type indicates with bounds will be constructed
/// @param bValues the bound values
///
/// @note will throw an exception if misconfigured, e.g. mismatching sType/bType
///       or mismatching bValues to bType
///
/// @return a shared surface instance
std::shared_ptr<Surface> createSurface(
    const Transform3& transform, Surface::SurfaceType sType,
    SurfaceBounds::BoundsType bType,
    const std::vector<ActsScalar>& bValues) noexcept(false);

/// @brief Create a bounds object from bounds type and bounds parameters
///
/// @param bType is the surface bounds type
/// @param bValues are the surface bounds values
///
/// @note will throw an exeption is mismatching bounds type with bValue length
///
/// @return a shared bounds instance
std::shared_ptr<SurfaceBounds> createBounds(
    SurfaceBounds::BoundsType bType,
    const std::vector<ActsScalar>& bValues) noexcept(false);

/// @brief Bounds creation from type and bound values
/// @tparam bounds_t the bounds type to be contructed
/// @param bValues the bound values
/// @return a shared object of the given bounds type
template <typename bounds_t>
std::shared_ptr<bounds_t> createBounds(const std::vector<ActsScalar>& bValues) {
  const size_t kValues = bounds_t::BoundValues::eSize;
  std::array<ActsScalar, kValues> bArray;
  std::copy_n(bValues.begin(), kValues, bArray.begin());
  return std::make_shared<bounds_t>(bArray);

}

}  // namespace SurfaceFactory

}  // namespace Acts
