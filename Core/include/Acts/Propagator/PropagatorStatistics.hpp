// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"

namespace Acts {

/// @struct PropagatorStatistics
///
/// @brief A struct to hold statistics of the propagator
struct PropagatorStatistics {
  /// Statistics of the stepper
  StepperStatistics stepping;
  /// Statistics of the navigator
  NavigatorStatistics navigation;
};

}  // namespace Acts
