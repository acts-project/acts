// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Hashing/HashingTraining.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/AngleHelpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <cmath>
#include <cstdint>

#include <annoy/annoylib.h>
#include <annoy/kissrandom.h>

namespace Acts {

HashingTraining::HashingTraining(const Config& cfg) : m_cfg(cfg) {
  if (m_cfg.f <= 0) {
    throw std::invalid_argument("Invalid f, f must be positive");
  }
  if (m_cfg.annoySeed <= 0) {
    throw std::invalid_argument(
        "Invalid Annoy random seed, Annoy random seed must be positive");
  }
}

AnnoyModel HashingTraining::execute(
    const SpacePointContainer2& spacePoints) const {
  const unsigned int annoySeed = m_cfg.annoySeed;
  const std::int32_t f = m_cfg.f;

  auto annoyModel = AnnoyModel(f);

  annoyModel.set_seed(annoySeed);

  for (const auto& spacePoint : spacePoints) {
    const float x = spacePoint.x() / Acts::UnitConstants::mm;
    const float y = spacePoint.y() / Acts::UnitConstants::mm;

    const float phi = std::atan2(y, x);

    std::array<float, 2> vec{};
    if (f >= 1) {
      vec[0] = phi;
    }
    if (f >= 2) {
      const float z = spacePoint.z() / Acts::UnitConstants::mm;
      const float r2 = hypotSquare(x, y);
      const float rho = std::sqrt(r2 + z * z);
      const float theta = std::acos(z / rho);
      const float eta = AngleHelpers::etaFromTheta(theta);
      vec[1] = eta;
    }

    annoyModel.add_item(spacePoint.index(), vec.data());
  }

  const int nTrees = 2 * f;
  annoyModel.build(nTrees);

  return annoyModel;
}

}  // namespace Acts
