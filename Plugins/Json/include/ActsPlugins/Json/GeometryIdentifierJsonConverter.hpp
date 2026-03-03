// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

namespace GeometryIdentifierJsonConverter {
/// Encode the geometry identifier
/// @param geoId is the geometry identifier that will be encoded
/// @param compact if true, the raw value is stored
inline nlohmann::json encodeIdentifier(const GeometryIdentifier& geoId,
                                       bool compact = false) {
  nlohmann::json encoded;
  if (compact) {
    // store the raw value directly
    encoded = geoId.value();
    return encoded;
  }

  // only store non-zero identifiers
  if (geoId.volume() != 0u) {
    encoded["volume"] = geoId.volume();
  }
  if (geoId.boundary() != 0u) {
    encoded["boundary"] = geoId.boundary();
  }
  if (geoId.layer() != 0u) {
    encoded["layer"] = geoId.layer();
  }
  if (geoId.approach() != 0u) {
    encoded["approach"] = geoId.approach();
  }
  if (geoId.sensitive() != 0u) {
    encoded["sensitive"] = geoId.sensitive();
  }
  if (geoId.extra() != 0u) {
    encoded["extra"] = geoId.extra();
  }
  return encoded;
}

/// @brief  Decode a geometry identifier from a json object
/// @param encoded is the json object that carries the encoded identifier
/// @return a valid geometry Identifier
inline GeometryIdentifier decodeIdentifier(const nlohmann::json& encoded) {
  return GeometryIdentifier()
      .withVolume(encoded.value("volume", GeometryIdentifier::Value{0u}))
      .withBoundary(encoded.value("boundary", GeometryIdentifier::Value{0u}))
      .withLayer(encoded.value("layer", GeometryIdentifier::Value{0u}))
      .withApproach(encoded.value("approach", GeometryIdentifier::Value{0u}))
      .withSensitive(encoded.value("sensitive", GeometryIdentifier::Value{0u}))
      .withExtra(encoded.value("extra", GeometryIdentifier::Value{0u}));
}

}  // namespace GeometryIdentifierJsonConverter

/// Write to JSON
/// @param j JSON object to write to
/// @param geoId GeometryIdentifier to write
void to_json(nlohmann::json& j, const GeometryIdentifier& geoId);

/// Read from JSON
/// @param j JSON object to read from
/// @param geoId  GeometryIdentifier to fill
void from_json(const nlohmann::json& j, GeometryIdentifier& geoId);

/// @}
}  // namespace Acts
