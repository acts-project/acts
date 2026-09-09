// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <stdexcept>

namespace ActsPlugins::detail {

/// @brief Convert a GeoModel placement into an ACTS transform
///
/// @c GeoTrf::Transform3D is an @c Eigen::Affine3d, while ACTS requires rigid
/// transforms. GeoModel placements are built from translations and rotations,
/// so this only has to reject the scaled or sheared ones.
///
/// @param transform The GeoModel transform to convert
/// @throws std::invalid_argument if the linear part is not orthogonal
/// @return The equivalent ACTS transform
inline Acts::Transform3 convertTransform(
    const Acts::AffineTransform3& transform) {
  if (!Acts::isOrthogonal(transform.linear())) {
    throw std::invalid_argument(
        "GeoModel placement is not a rigid transformation");
  }
  return Acts::makeTransform3(transform.linear(), transform.translation());
}

}  // namespace ActsPlugins::detail
