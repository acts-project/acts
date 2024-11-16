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

/// @struct NavigatorStatistics
///
/// @brief A struct to hold statistics of the navigator
struct NavigatorStatistics {
  /// Number of renavigation attempts
  std::size_t nRenavigations = 0;

  /// Number of volume switches
  std::size_t nVolumeSwitches = 0;
};

}  // namespace Acts
