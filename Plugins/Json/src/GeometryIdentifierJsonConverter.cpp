// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"

void Acts::to_json(nlohmann::json& j, const GeometryIdentifier& geoId) {
  j = GeometryIdentifierJsonConverter::encodeIdentifier(geoId);
}

void Acts::from_json(const nlohmann::json& j, GeometryIdentifier& geoId) {
  // Check if the JSON object iself is a GeometryIdentifier::ValueType
  if (!j.is_object()) {
    auto value = j.get<GeometryIdentifier::Value>();
    geoId = GeometryIdentifier(value);
    return;
  }

  geoId = GeometryIdentifierJsonConverter::decodeIdentifier(j);
}
