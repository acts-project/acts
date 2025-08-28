// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.

namespace Acts {
class BinningData;

/// Convert BinningData to JSON
/// @param j Destination JSON object
/// @param bd Source BinningData to convert
void to_json(nlohmann::json& j, const BinningData& bd);

/// Convert JSON to BinningData
/// @param j Source JSON object
/// @param bd Destination BinningData to populate
void from_json(const nlohmann::json& j, BinningData& bd);

/// Convert BinUtility to JSON
/// @param j Destination JSON object
/// @param bu Source BinUtility to convert
void to_json(nlohmann::json& j, const BinUtility& bu);

/// Convert JSON to BinUtility
/// @param j Source JSON object
/// @param bu Destination BinUtility to populate
void from_json(const nlohmann::json& j, BinUtility& bu);

/// Convert Range1D to JSON
/// @param j Destination JSON object
/// @param r Source Range1D to convert
template <typename Type>
void to_json(nlohmann::json& j, const Range1D<Type>& r) {
  j["min"] = r.min();
  j["max"] = r.max();
}

/// Convert JSON to Range1D
/// @param j Source JSON object
/// @param r Destination Range1D to populate
template <typename Type>
void from_json(const nlohmann::json& j, Range1D<Type>& r) {
  r.setMin(static_cast<Type>(j["min"]));
  r.setMax(static_cast<Type>(j["max"]));
}

NLOHMANN_JSON_SERIALIZE_ENUM(AxisDirection,
                             {{AxisDirection::AxisX, "AxisX"},
                              {AxisDirection::AxisY, "AxisY"},
                              {AxisDirection::AxisZ, "AxisZ"},
                              {AxisDirection::AxisR, "AxisR"},
                              {AxisDirection::AxisPhi, "AxisPhi"},
                              {AxisDirection::AxisRPhi, "AxisRPhi"},
                              {AxisDirection::AxisTheta, "AxisTheta"},
                              {AxisDirection::AxisEta, "AxisEta"},
                              {AxisDirection::AxisMag, "AxisMag"}})

}  // namespace Acts
