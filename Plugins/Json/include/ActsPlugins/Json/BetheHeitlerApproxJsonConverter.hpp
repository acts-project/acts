// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"

#include <string>

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

void to_json(nlohmann::json& j,
             const PolynomialBetheHeitlerApprox::PolyData& data);

void from_json(const nlohmann::json& j,
               PolynomialBetheHeitlerApprox::PolyData& data);

void to_json(nlohmann::json& j,
             const PolynomialBetheHeitlerApprox::RangeData& data);

void from_json(const nlohmann::json& j,
               PolynomialBetheHeitlerApprox::RangeData& data);

/// Load a Bethe-Heitler approximation from a JSON file.
///
/// @param filepath Path to the JSON file
/// @param clampToRange Whether to clamp x/x0 values to the valid range
/// @param noChangeLimit Limit below which no change is applied
/// @param singleGaussianLimit Limit below which a single Gaussian is used
/// @return The loaded PolynomialBetheHeitlerApprox
PolynomialBetheHeitlerApprox loadBetheHeitlerApproxFromJson(
    const std::string& filepath, bool clampToRange = false,
    double noChangeLimit = 0.0001, double singleGaussianLimit = 0.002);

/// @}

}  // namespace Acts
