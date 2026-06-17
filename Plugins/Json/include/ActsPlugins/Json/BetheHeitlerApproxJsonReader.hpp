// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"

#include <memory>
#include <string>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// Load a Bethe-Heitler approximation from a JSON file.
///
/// @param filepath Path to the JSON file
/// @param clampToRange Whether to clamp x/x0 values to the valid range
/// @param noChangeLimit Limit below which no change is applied (default: 0.0001)
/// @param singleGaussianLimit Limit below which a single Gaussian is used (default: 0.002)
/// @return Shared pointer to the loaded PolynomialBetheHeitlerApprox
///
/// The JSON file should contain:
/// - "version": string (optional)
/// - "ranges": array of range objects with:
///   - "low_x0": lower limit of x/x0 range
///   - "high_x0": upper limit of x/x0 range
///   - "transform": whether to apply logistic transform (default: true)
///   - "components": array of component objects with "weight_coeffs",
///   "mean_coeffs", "var_coeffs"
std::shared_ptr<const PolynomialBetheHeitlerApprox>
loadBetheHeitlerApproxFromJson(const std::string& filepath,
                               bool clampToRange = false,
                               double noChangeLimit = 0.0001,
                               double singleGaussianLimit = 0.002);

/// @}

}  // namespace Acts
