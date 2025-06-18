// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <numbers>
#include <random>

namespace ActsExamples {

/// @brief Generator interface for a vertex position
struct PrimaryVertexPositionGenerator {
  /// @brief Virtual destructor required
  virtual ~PrimaryVertexPositionGenerator() = default;
  /// @brief Generate a vertex position
  ///
  /// @param rng Shared random number generator instance
  /// @return Acts::Vector4 The vertex position
  virtual Acts::Vector4 operator()(RandomEngine& rng) const = 0;
};

struct FixedPrimaryVertexPositionGenerator
    : public PrimaryVertexPositionGenerator {
  /// The fixed vertex position and time.
  Acts::Vector4 fixed = Acts::Vector4::Zero();

  Acts::Vector4 operator()(RandomEngine& /*rng*/) const override {
    return fixed;
  }
};

struct GaussianPrimaryVertexPositionGenerator
    : public PrimaryVertexPositionGenerator {
  // standard deviation comes first, since it is more likely to be modified
  /// Vertex position and time standard deviation.
  Acts::Vector4 stddev = {0.0, 0.0, 0.0, 0.0};
  /// Mean vertex position and time.
  Acts::Vector4 mean = {0.0, 0.0, 0.0, 0.0};

  Acts::Vector4 operator()(RandomEngine& rng) const override {
    auto normal = std::normal_distribution<double>(0.0, 1.0);
    Acts::Vector4 rndNormal = {
        normal(rng),
        normal(rng),
        normal(rng),
        normal(rng),
    };
    return mean + rndNormal.cwiseProduct(stddev);
  }
};

//
struct GaussianDisplacedVertexPositionGenerator
    : public PrimaryVertexPositionGenerator {
  double rMean = 0;
  double rStdDev = 1;
  double zMean = 0;
  double zStdDev = 1;
  double tMean = 0;
  double tStdDev = 1;

  Acts::Vector4 operator()(RandomEngine& rng) const override {
    double min_value = -std::numbers::pi;
    double max_value = std::numbers::pi;

    std::uniform_real_distribution uniform(min_value, max_value);

    std::normal_distribution rDist(rMean, rStdDev);
    std::normal_distribution zDist(zMean, zStdDev);
    std::normal_distribution tDist(tMean, tStdDev);

    // Generate random values from normal distributions
    double r = rDist(rng);
    double phi = uniform(rng);  // Random angle in radians
    double z = zDist(rng);
    double t = tDist(rng);

    // Convert cylindrical coordinates to Cartesian coordinates
    double x = r * std::cos(phi);
    double y = r * std::sin(phi);

    return Acts::Vector4(x, y, z, t);
  }
};
}  // namespace ActsExamples
