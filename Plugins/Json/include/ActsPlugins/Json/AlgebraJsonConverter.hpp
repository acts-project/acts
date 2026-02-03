// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// Convert Transform3 to JSON
/// @param j Destination JSON object
/// @param t Source Transform3 to convert
void to_json(nlohmann::json& j, const Transform3& t);

/// Convert JSON to Transform3
/// @param j Source JSON object
/// @param t Destination Transform3 to populate
void from_json(const nlohmann::json& j, Transform3& t);

namespace Transform3JsonConverter {

/// @brief The options for the transform converter
struct Options {
  /// Write the identity transform explicitly
  bool writeIdentity = false;
  /// Apply a transpose to flip column/row order
  bool transpose = false;
};

/// @brief The Transform converter to json
///
/// @param t the transform to be converted
/// @param options transformation options
///
/// @return a json object representing the transform
nlohmann::json toJson(const Transform3& t, const Options& options = {});

/// @brief The Transform converter from json
///
/// @param jTransform the transform json transformation
///
/// @return a transform object
Transform3 fromJson(const nlohmann::json& jTransform);

}  // namespace Transform3JsonConverter

/// @}
}  // namespace Acts
