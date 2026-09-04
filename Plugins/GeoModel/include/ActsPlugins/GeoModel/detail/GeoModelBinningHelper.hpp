// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisDefinitions.hpp"

#include <stdexcept>
#include <string>

namespace ActsPlugins::detail::GeoModelBinningHelper {

/// @brief Helper to transform binning string to AxisDirection enum
///
/// @param binning the binning string
inline Acts::AxisDirection toAxisDirection(const std::string& binning) {
  using enum Acts::AxisDirection;
  if (binning == "x") {
    return AxisX;
  } else if (binning == "y") {
    return AxisY;
  } else if (binning == "z") {
    return AxisZ;
  } else if (binning == "r") {
    return AxisR;
  } else if (binning == "phi") {
    return AxisPhi;
  }
  throw std::invalid_argument("GeoModelBinningHelper: Unknown binning value '" +
                              binning + "'");
}

}  // namespace ActsPlugins::detail::GeoModelBinningHelper
