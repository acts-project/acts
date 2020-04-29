// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2018-03-13
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include <array>
#include <random>

#include "ACTFW/Framework/RandomNumbers.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace FW {

struct FixedVertexGenerator {
  /// The fixed vertex position and time.
  Acts::ActsVector<double, 4> fixed = {0.0, 0.0, 0.0, 0.0};

  Acts::ActsVector<double, 4> operator()(RandomEngine& /* unused */) const {
    return fixed;
  }
};

struct GaussianVertexGenerator {
  // standard deviation comes first, since it is more likely to be modified
  /// Vertex position and time standard deviation.
  std::array<double, 4> stddev = {0.0, 0.0, 0.0, 0.0};
  /// Mean vertex position and time.
  std::array<double, 4> mean = {0.0, 0.0, 0.0, 0.0};

  Acts::ActsVector<double, 4> operator()(RandomEngine& rng) const {
    auto x = std::normal_distribution<double>(mean[0], stddev[0])(rng);
    auto y = std::normal_distribution<double>(mean[1], stddev[1])(rng);
    auto z = std::normal_distribution<double>(mean[2], stddev[2])(rng);
    auto t = std::normal_distribution<double>(mean[3], stddev[3])(rng);
    return {x, y, z, t};
  }
};

}  // namespace FW
