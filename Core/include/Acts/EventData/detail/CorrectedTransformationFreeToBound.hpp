// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"

namespace Acts {
namespace detail {

/// @brief Parameters sampling based on covariance matrix sqrt root in UKF: https://doi.org/10.1117/12.280797
///
///@todo make the parameters configurable
struct CorrectedFreeToBoundTransformer {
  /// The parameter to tune the weight (0 < alpha <=1)
  ActsScalar alpha = 0.1;
  ActsScalar beta = 2;

  // The incident angle cutoff for correction
  ActsScalar incidentAngleCutoff = 0.1;

  /// Get the non-linearity corrected bound parameters and its covariance
  std::optional<std::tuple<BoundVector, BoundSymMatrix, BoundToFreeMatrix>>
  operator()(const FreeVector& freeParams, const FreeSymMatrix& freeCovariance,
             const Surface& surface, const GeometryContext& geoContext,
             NavigationDirection navDir = NavigationDirection::Forward);
};

}  // namespace detail
}  // namespace Acts
