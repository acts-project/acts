// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

namespace BetheHeitlerApproxJsonConverter {

/// Single x/x0 range data for one component
struct RangeData {
  double lowX0 = 0;
  double highX0 = 0;
  std::vector<double> weightCoeffs;
  std::vector<double> meanCoeffs;
  std::vector<double> varCoeffs;
};

/// Convert PolyData (component coefficients) to JSON
/// @param j Destination JSON object
/// @param data Source PolyData to convert
void to_json(nlohmann::json& j,
             const PolynomialBetheHeitlerApprox::PolyData& data);

/// Convert JSON to PolyData (component coefficients)
/// @param j Source JSON object
/// @param data Destination PolyData to populate
void from_json(const nlohmann::json& j,
               PolynomialBetheHeitlerApprox::PolyData& data);

/// Convert Data (vector of PolyData) to JSON for a single x0 range
/// @param j Destination JSON object
/// @param data Source Data to convert
void to_json(nlohmann::json& j, const PolynomialBetheHeitlerApprox::Data& data);

/// Convert JSON to Data (vector of PolyData) for a single x0 range
/// @param j Source JSON object
/// @param data Destination Data to populate
void from_json(const nlohmann::json& j,
               PolynomialBetheHeitlerApprox::Data& data);

/// Convert RangeData to JSON
/// @param j Destination JSON object
/// @param data Source RangeData to convert
void to_json(nlohmann::json& j, const RangeData& data);

/// Convert JSON to RangeData
/// @param j Source JSON object
/// @param data Destination RangeData to populate
void from_json(const nlohmann::json& j, RangeData& data);

}  // namespace BetheHeitlerApproxJsonConverter

/// @}

}  // namespace Acts
