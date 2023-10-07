// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Range1D.hpp"

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.

namespace Acts {

void to_json(nlohmann::json& j, const BinningData& bd);

void from_json(const nlohmann::json& j, BinningData& bd);

void to_json(nlohmann::json& j, const BinUtility& bu);

void from_json(const nlohmann::json& j, BinUtility& bu);

template <typename Type>
void to_json(nlohmann::json& j, const Range1D<Type>& r) {
  j["min"] = r.min();
  j["max"] = r.max();
}

template <typename Type>
void from_json(const nlohmann::json& j, Range1D<Type>& r) {
  r.setMin(static_cast<Type>(j["min"]));
  r.setMax(static_cast<Type>(j["max"]));
}

}  // namespace Acts
