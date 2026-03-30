// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/ProtoAxis.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

/// @namespace Acts::ProtoAxisJsonConverter
/// @ingroup json_plugin

namespace Acts::ProtoAxisJsonConverter {

/// @addtogroup json_plugin
/// @{

/// Write the ProtoAxis to a json object
///
/// @param pa the proto axis to be written out
/// @return JSON object representing the proto axis
nlohmann::json toJson(const ProtoAxis& pa);

/// Create a ProtoAxis from a json object
///
/// @param j the json object to be read from
/// @return ProtoAxis created from the JSON object
Acts::ProtoAxis fromJson(const nlohmann::json& j);

/// @}

}  // namespace Acts::ProtoAxisJsonConverter
