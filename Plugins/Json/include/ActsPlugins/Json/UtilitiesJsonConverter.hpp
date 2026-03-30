// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{
class BinningData;

/// Convert BinningData to JSON
/// @param j JSON object to write to
/// @param bd BinningData to convert
void to_json(nlohmann::json& j, const BinningData& bd);

/// Convert JSON to BinningData
/// @param j The JSON object to convert from
/// @param bd The BinningData object to populate
void from_json(const nlohmann::json& j, BinningData& bd);

/// Convert BinUtility to JSON
/// @param j JSON object to write to
/// @param bu BinUtility to convert
void to_json(nlohmann::json& j, const BinUtility& bu);

/// Convert JSON to BinUtility
/// @param j JSON object to convert from
/// @param bu BinUtility to populate
void from_json(const nlohmann::json& j, BinUtility& bu);

/// Convert Range1D to JSON
/// @param j JSON object to write to
/// @param r Range1D to convert
template <typename Type>
void to_json(nlohmann::json& j, const Range1D<Type>& r) {
  j["min"] = r.min();
  j["max"] = r.max();
}

/// Convert JSON to Range1D
/// @param j JSON object to convert from
/// @param r Range1D to populate
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

/// @}
}  // namespace Acts
