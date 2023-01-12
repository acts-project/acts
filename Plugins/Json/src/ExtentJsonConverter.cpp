// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/ExtentJsonConverter.hpp"

#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Enumerate.hpp"

void Acts::to_json(nlohmann::json& j, const Acts::Extent& e) {
  const auto bValueNames = binningValueNames();
  const auto& xrange = e.range();
  for (auto [ib, ibv] : enumerate(s_binningValues)) {
    if (e.constrains(ibv)) {
      const auto& r = xrange[ib];
      j[bValueNames[ib]] = r;
    }
  }
}

void Acts::from_json(const nlohmann::json& j, Acts::Extent& e) {
  const auto bValueNames = binningValueNames();
  for (auto [ib, bvn] : enumerate(bValueNames)) {
    if (j.find(bvn) != j.end()) {
      e.set(static_cast<BinningValue>(ib), j[bvn]["min"], j[bvn]["max"]);
    }
  }
}
