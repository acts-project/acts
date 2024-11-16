// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>

namespace Acts {

/// @struct StepperStatistics
///
/// @brief A struct to hold statistics of the stepper
struct StepperStatistics {
  /// Number of attempted steps
  std::size_t nAttemptedSteps = 0;
  /// Number of rejected steps
  std::size_t nRejectedSteps = 0;
  /// Number of successful steps
  std::size_t nSuccessfulSteps = 0;
  /// Number of steps that were reversed
  std::size_t nReverseSteps = 0;

  /// Signed sum of the step lengths
  double pathLength = 0;
  /// Unsigned sum of the step lengths
  double absolutePathLength = 0;
};

}  // namespace Acts
