// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"

#include <random>

namespace ActsExamples {

struct FixedPrimaryVertexPositionGenerator
    : public EventGenerator::PrimaryVertexPositionGenerator {
  /// The fixed vertex position and time.
  Acts::Vector4 fixed = Acts::Vector4::Zero();

  Acts::Vector4 operator()(RandomEngine& /*rng*/) const override {
    return fixed;
  }
};

struct GaussianPrimaryVertexPositionGenerator
    : public EventGenerator::PrimaryVertexPositionGenerator {
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
    : public EventGenerator::PrimaryVertexPositionGenerator {
  Acts::Vector4 stddev = {0.0, 0.0, 0.0, 0.0};
  Acts::Vector4 mean = {0.0, 0.0, 0.0, 0.0};

  Acts::Vector4 operator()(RandomEngine& rng) const override {
    double min_value = -2 * M_PI;  // -π
    double max_value = 2 * M_PI;   // π

    // Create a uniform real distribution object with range [-π, π)
    std::uniform_real_distribution<> uniform(min_value, max_value);

    // Create a normal distribution object with mean 0.0 and standard
    // deviation 1.0
    auto normal = std::normal_distribution<double>(0.0, 1.0);
    Acts::Vector4 rnd = {normal(rng), uniform(rng), normal(rng), normal(rng)};

    // Compute the cylindrical coordinates with deviations applied
    Acts::Vector4 cylindrical = mean + rnd.cwiseProduct(stddev);

    // Extract coordinates from the vector
    double r = cylindrical(0);
    double phi = cylindrical(1);
    double z = cylindrical(2);
    double t = cylindrical(3);
    double phi_degrees = phi * (180.0 / M_PI);

    // Convert cylindrical coordinates to Cartesian coordinates
    double x = r * std::cos(phi_degrees);
    double y = r * std::sin(phi_degrees);

    return Acts::Vector4(x, y, z, t);
  }
};

}  // namespace ActsExamples
