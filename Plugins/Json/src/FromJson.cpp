// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/lib/json.hpp"
#include "Acts/Utilities/VariantData.hpp"

#include "Acts/Plugins/Json/FromJson.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

using json = nlohmann::json;

Acts::variant_data
Acts::from_json(const json& node)
{
  throw_assert(!node.is_null(), "Null json values cannot be deserialized");

  if (node.is_array()) {
    variant_vector vec;
    for (const auto& item : node) {
      vec.push_back(from_json(item));
    }

    return vec;
  } else if (node.is_object()) {

    variant_map map;
    for (json::const_iterator it = node.begin(); it != node.end(); ++it) {
      map[it.key()] = from_json(it.value());
    }

    return map;
  } else if (node.is_boolean()) {
    bool value = node;
    return value;
  } else if (node.is_string()) {
    std::string value = node;
    return value;
  } else if (node.is_number_float()) {
    double value = node;
    return value;
  } else /*if(node.is_number_integer())*/ {
    int value = node;
    return value;
  }
}
