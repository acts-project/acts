// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Propagator/PropagatorStatistics.hpp"

#include <optional>

namespace Acts {

/// @brief Simple class holding result of propagation call
///
/// The actor state stays in the propagator state. A caller that needs it runs
/// the propagation on its own state and reads it there.
///
/// @tparam parameters_t Type of final track parameters
template <typename parameters_t>
struct PropagatorResult {
  /// Final track parameters
  std::optional<parameters_t> endParameters = std::nullopt;

  /// Full transport jacobian
  std::optional<BoundMatrix> transportJacobian = std::nullopt;

  /// Number of propagation steps that were carried out
  std::size_t steps = 0;

  /// Signed distance over which the parameters were propagated
  double pathLength = 0.;

  /// Propagator statistics
  PropagatorStatistics statistics;
};

}  // namespace Acts
