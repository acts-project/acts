// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

namespace ActsPlugins::detail {

/// @brief Convert a GeoModel placement into an ACTS transform
///
/// GeoModel stores placements as general affine transforms
/// (@c GeoTrf::Transform3D is an @c Eigen::Affine3d), while ACTS requires them
/// to be rigid. Placements in a GeoModel tree are built from translations and
/// rotations, so the orthogonality of the linear part is only asserted here
/// rather than enforced.
///
/// @param transform The GeoModel transform to convert
/// @return The equivalent ACTS transform
inline Acts::Transform3 convertTransform(
    const Acts::AffineTransform3& transform) {
  return Acts::makeTransform3(transform.linear(), transform.translation());
}

}  // namespace ActsPlugins::detail
