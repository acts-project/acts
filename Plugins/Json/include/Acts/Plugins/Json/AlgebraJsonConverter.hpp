// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"

#include <unordered_map>

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.
namespace Acts {

void to_json(nlohmann::json& j, const Transform3& t);

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

/// Read the Transform converter from json
///
/// @param jTransform the transform json transformation
///
/// @return a transform object
Transform3 fromJson(const nlohmann::json& jTransform);

}  // namespace Transform3JsonConverter

namespace IdentifiedTransform3JsonConverter {

/// @brief The options for the transform converter
struct Options {
  /// How the transform is written
  Transform3JsonConverter::Options transformOptions;
  /// Write the geometry identifier only as a single value
  bool compressIdentifier = false;
};

using Collection = std::unordered_map<GeometryIdentifier, Transform3>;

/// The identified transform container json I/O
///
/// @param identifiedTransforms is the map with geometry id and transforms
/// @param options transformation options
///
/// @return a json object representing this map
nlohmann::json toJson(const Collection& identifiedTransforms,
                      const Options& options = {});

/// Read the Identified transforms from json
///
/// @param jIdentifiedTransform the json representing the container
///
/// @return an unordered map with identified transforms
Collection fromJson(const nlohmann::json& jIdentifiedTransforms);

}  // namespace IdentifiedTransform3JsonConverter
}  // namespace Acts
