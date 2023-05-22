// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <array>
#include <vector>

namespace Acts {
namespace Experimental {
namespace detail {
/// @brief  This file contains helper methods to build common support structures
/// such as support cylinders or discs.
///
/// It allows to model those as Disc/CylinderSurface objects, but also - if
/// configured such - as approximations built from palanr surfaces
namespace SupportHelper {

/// @brief Helper method to build cylindrical support structure
///
/// @param transform the transform where this cylinder should be located
/// @param bounds the boundary values: r, halfZ, halfPhi, avgPhi, bevell parameters
/// @param splits the number of surfaces through which the surface is approximated (1u ... cylinder)
///
/// @return a vector of surfaces that represent this support
std::vector<std::shared_ptr<Surface>> cylindricalSupport(
    const Transform3& transform, const std::array<ActsScalar, 6u>& bounds,
    unsigned int splits = 1u);

/// @brief Helper method to build disc support structure
///
/// @param transform the transform where this disc should be located
/// @param bounds the boundary values: minR, maxR, halfPhi, avgPhi
/// @param splits the number of surfaces through which the surface is approximated (1u ... disc)
///
/// @return a vector of surfaces that represent this support
std::vector<std::shared_ptr<Surface>> discSupport(
    const Transform3& transform, const std::array<ActsScalar, 4u>& bounds,
    unsigned int splits = 1u);

/// Add support to already existing surfaces
///
/// @param layerSurfaces [in, out] the surfaces to which those are added
/// @param assignToAll [in, out] indices that are assigned to all bins in the indexing
/// @param layerExtent the externally provided layer Extent
/// @param layerRepresentation the shape of the layer
/// @param layerSupportValues the support structue in numbers
/// @param layerTransform is an optional value of the layer transform
/// @param supportSplits the number of splits if splitting is configured
///
/// @note this modifies the layerSurfaces and toAllIndices
void addSupport(std::vector<std::shared_ptr<Surface>>& layerSurfaces,
                std::vector<size_t>& assignToAll, const Extent& layerExtent,
                Surface::SurfaceType layerRepresentation,
                const std::array<ActsScalar, 5u>& layerSupportValues,
                std::optional<Transform3> layerTransform = std::nullopt,
                unsigned int supportSplits = 1u);

}  // namespace SupportHelper
}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
