// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/ExtentJsonConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <array>
#include <iterator>
#include <vector>

void Acts::to_json(nlohmann::json& j, const Acts::Extent& e) {
  const auto& bValueNames = binningValueNames();

  {
    nlohmann::json jrange;
    const auto& xrange = e.range();
    for (auto [ib, ibv] : enumerate(s_binningValues)) {
      if (e.constrains(ibv)) {
        jrange[bValueNames[ib]] = xrange[ib];
      }
    }
    j["range"] = jrange;
  }

  {
    nlohmann::json jenvelope;
    const auto& envelope = e.envelope();
    for (auto [ib, ibv] : enumerate(s_binningValues)) {
      if (envelope[ibv] != zeroEnvelope) {
        jenvelope[bValueNames[ib]] =
            Range1D<ActsScalar>(envelope[ib][0], envelope[ib][1]);
      }
    }
    if (!jenvelope.empty()) {
      j["envelope"] = jenvelope;
    }
  }
}

void Acts::from_json(const nlohmann::json& j, Acts::Extent& e) {
  const auto& bValueNames = binningValueNames();

  {
    const auto& jrange = j["range"];
    for (auto [ib, bvn] : enumerate(bValueNames)) {
      if (jrange.contains(bvn)) {
        e.set(static_cast<BinningValue>(ib), jrange[bvn]["min"],
              jrange[bvn]["max"]);
      }
    }
  }

  if (j.contains("envelope")) {
    const auto& jenvelope = j["envelope"];
    ExtentEnvelope envelope;
    for (auto [ib, bvn] : enumerate(bValueNames)) {
      if (jenvelope.find(bvn) != jenvelope.end()) {
        envelope[ib] = {jenvelope[bvn]["min"], jenvelope[bvn]["max"]};
      }
    }
    e.setEnvelope(envelope);
  }
}
