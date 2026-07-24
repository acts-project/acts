// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/AxisFactory.hpp"
#include "Acts/Utilities/MultiAxisFactory.hpp"

#include <memory>
#include <vector>

namespace Acts {

class Surface;

/// @brief Get the canonical local axis directions of a surface, in order
///
/// - Cylinder: (rphi, z)
/// - Disc with radial bounds: (r, phi)
/// - Plane with rectangle or trapezoid bounds: (x, y)
///
/// @param surface the surface to inspect
/// @throws std::invalid_argument for unsupported surface types
/// @return the axis directions, in canonical order
std::vector<AxisDirection> surfaceAxisDirections(const Surface& surface);

/// @brief Get the range and boundary type of a surface along a given axis
/// direction
///
/// The range is derived from the surface bounds in local coordinates.
/// Azimuthal directions on surfaces that cover the full azimuth resolve to a
/// closed boundary type; all other combinations resolve to a bound boundary
/// type.
///
/// @param surface the surface to inspect
/// @param aDir the axis direction to resolve
/// @throws std::invalid_argument for unsupported surface type and direction
///         combinations
/// @return the range and boundary type along the given direction
AxisResolution surfaceAxisResolution(const Surface& surface,
                                     AxisDirection aDir);

/// @brief Resolve a multi-dimensional binning description against a surface
///
/// Each axis description is matched to one of the canonical axis directions
/// of the surface (see @c surfaceAxisDirections ):
///
/// - if no description carries a direction, they are matched positionally,
///   i.e. slot i binds to the i-th canonical direction
/// - if all descriptions carry directions, each has to be one of the
///   canonical directions and no direction may appear twice; descriptions
///   given in a different order are re-ordered to match the surface
/// - a mixture of descriptions with and without directions is rejected
///
/// The number of descriptions may be smaller than the number of canonical
/// directions (e.g. a 1D binning on a cylinder). The resolved axes are
/// returned in canonical direction order, each carrying its direction.
///
/// @param multiAxisFactory the binning description to resolve
/// @param surface the surface to resolve against
/// @throws std::invalid_argument if the directions cannot be matched or the
///         surface is unsupported
/// @throws std::domain_error if any description is fully specified instead
///         of deferred
/// @return the resolved axes, in canonical direction order
std::vector<std::unique_ptr<IAxis>> resolveAxes(
    const MultiAxisFactory& multiAxisFactory, const Surface& surface);

}  // namespace Acts
