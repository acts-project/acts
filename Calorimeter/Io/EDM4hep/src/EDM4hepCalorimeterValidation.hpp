// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "ActsCalorimeter/EventData/CalorimeterHit.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <edm4hep/CalorimeterHitCollection.h>

namespace ActsExamples::detail {
inline float calorimeterFloat(double value) {
  if (!std::isfinite(value) ||
      std::abs(value) > std::numeric_limits<float>::max()) {
    throw std::overflow_error("Calorimeter value exceeds EDM4hep float range");
  }
  return static_cast<float>(value);
}
inline float calorimeterEnergy(double value) {
  const auto result = calorimeterFloat(value / Acts::UnitConstants::GeV);
  if (result <= 0) {
    throw std::overflow_error(
        "Positive calorimeter energy underflows EDM4hep float range");
  }
  return result;
}
inline bool calorimeterClose(double first, double second, double scale) {
  return std::isfinite(first) && std::isfinite(second) &&
         std::abs(first - second) <=
             64 * std::numeric_limits<double>::epsilon() * scale;
}
inline void checkCalorimeterHitMapping(
    const ActsCalorimeter::CalorimeterHitContainer& native,
    const edm4hep::CalorimeterHitCollection& persistent) {
  if (native.size() != persistent.size()) {
    throw std::invalid_argument(
        "Native/persistent hit collection sizes disagree");
  }
  for (std::size_t index = 0; index < native.size(); ++index) {
    const auto& hit = native[index];
    const auto output = persistent[index];
    if (!std::isfinite(hit.energy) || hit.energy <= 0 ||
        !hit.position.allFinite() || !std::isfinite(hit.time) ||
        !output.isAvailable() || output.getCellID() != hit.cellId ||
        output.getEnergy() != calorimeterEnergy(hit.energy) ||
        output.getTime() !=
            calorimeterFloat(hit.time / Acts::UnitConstants::ns) ||
        output.getPosition().x !=
            calorimeterFloat(hit.position.x() / Acts::UnitConstants::mm) ||
        output.getPosition().y !=
            calorimeterFloat(hit.position.y() / Acts::UnitConstants::mm) ||
        output.getPosition().z !=
            calorimeterFloat(hit.position.z() / Acts::UnitConstants::mm)) {
      throw std::invalid_argument(
          "Native/persistent hit mapping is inconsistent");
    }
  }
}
}  // namespace ActsExamples::detail
