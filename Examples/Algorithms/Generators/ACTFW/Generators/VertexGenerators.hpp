// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2018-03-13
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include <random>

#include "ACTFW/Framework/RandomNumbers.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace FW {

struct FixedVertexGenerator {
  /// The fixed vertex position and time.
  Acts::Vector4D fixed = Acts::Vector4D::Zero();

  Acts::Vector4D operator()(RandomEngine& /* unused */) const { return fixed; }
};

struct GaussianVertexGenerator {
  // standard deviation comes first, since it is more likely to be modified
  /// Vertex position and time standard deviation.
  Acts::Vector4D stddev = {0.0, 0.0, 0.0, 0.0};
  /// Mean vertex position and time.
  Acts::Vector4D mean = {0.0, 0.0, 0.0, 0.0};

  Acts::Vector4D operator()(RandomEngine& rng) const {
    auto normal = std::normal_distribution<double>(0.0, 1.0);
    Acts::Vector4D rndNormal = {
        normal(rng),
        normal(rng),
        normal(rng),
        normal(rng),
    };
    return mean + rndNormal.cwiseProduct(stddev);
  }
};

}  // namespace FW
