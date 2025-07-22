// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <functional>
#include <random>

namespace ActsExamples::AlignmentGenerator {

/// @note Non-contextual random number generators
///
/// Random numbers in the alignment generators can be context-free
/// as their call sequence is controlled  by the IOV and not altered
/// by eventually changing random number seeds.

struct UniformRandom {
  /// The random number generator to be used for normal distribution
  /// @param seed The random seed to be used
  /// @param min The minimum value of the uniform distribution
  /// @param max The maximum value of the uniform distribution
  UniformRandom(RandomSeed seed, double min, double max)
      : randomEngine(seed), randomNumbers(min, max) {}

  /// @brief The call operation generating a random number
  double operator()() { return randomNumbers(randomEngine); }

  RandomEngine randomEngine;  ///< The random number generator
  std::uniform_real_distribution<double> randomNumbers;  ///< The distribution
};

struct GaussRandom {
  /// The random number generator to be used for Gaussian distribution
  /// @param seed The random seed to be used
  /// @param mean The mean value of the Gaussian distribution
  /// @param sigma The standard deviation of the Gaussian distribution
  GaussRandom(RandomSeed seed, double mean, double sigma)
      : randomEngine(seed), randomNumbers(mean, sigma) {}

  /// @brief The call operation generating a random number
  double operator()() { return randomNumbers(randomEngine); }

  RandomEngine randomEngine;  ///< The random number generator
  std::normal_distribution<double> randomNumbers;  ///< The distribution
};

struct ClippedGaussRandom {
  /// The random number generator to be usedfor gaussian with clipping
  /// @param seed The random seed to be used
  /// @param mean The mean value of the Gaussian distribution
  /// @param sigma The standard deviation of the Gaussian distribution
  /// @param minv The minimum value of the truncated Gaussian distribution
  /// @param maxv The maximum value of the truncated Gaussian distribution
  ClippedGaussRandom(RandomSeed seed, double mean, double sigma, double minv,
                     double maxv)
      : randomEngine(seed), randomNumbers(mean, sigma), min(minv), max(maxv) {}

  /// @brief The call operation generating a random number
  double operator()() {
    auto value = randomNumbers(randomEngine);
    return std::max(min, std::min(max, value));
  }

  RandomEngine randomEngine;  ///< The random number generator
  std::normal_distribution<double> randomNumbers;  ///< The distribution
  double min;  ///< The minimum value of the truncated Gaussian distribution
  double max;  ///< The maximum value of the truncated Gaussian distribution
};

/// This generator does nothing
struct Nominal {
  /// @brief The call operation for the nominal alignment
  void operator()(Acts::Transform3* /*transform*/) const {
    // No operation, this is a nominal alignment generator}
  }
};

/// This generator applies a constant global shift to the transform
struct GlobalShift {
  // The shift to be applied
  Acts::Vector3 shift = Acts::Vector3::UnitZ();

  std::function<double()> randomize = nullptr;  ///< Optional scale factor

  /// @brief The call operation applying the global shift
  void operator()(Acts::Transform3* transform) const {
    double scale =
        randomize != nullptr ? randomize() : 1.0;  ///< Default scale is 1.0
    transform->translation() += scale * shift;
  }
};

/// This generator applies a constant global rotation (i.e. around 0.,0.,0)
/// of value `angle` around the `axis`
struct GlobalRotation {
  /// The axis around which the rotation is applied
  Acts::Vector3 axis =
      Acts::Vector3::UnitZ();  ///< The rotation axis, default is Z-axis
  double angle = 0.0;          ///< The rotation angle in radians

  std::function<double()> randomize = nullptr;  ///< Optional scale factor

  /// @brief The call operation applying the global rotation
  /// @param transform The transform to be rotated
  void operator()(Acts::Transform3* transform) const {
    double angleS = randomize != nullptr ? randomize() * angle : angle;
    (*transform) = Acts::AngleAxis3(angleS, axis) * (*transform);
  }
};

/// This generator applies a rotation in the local frame around a given
/// local axis. I.e. if the z-axis is chosen, the module will be rotated
/// by the angle around this axis.
struct LocalRotation {
  /// The axis around which the rotation is applied
  Acts::Vector3 axis =
      Acts::Vector3::UnitZ();  ///< The rotation axis, default is Z-axis
  double angle = 0.0;          ///< The rotation angle in radians

  std::function<double()> randomize = nullptr;  ///< Optional scale factor

  /// @brief The call operation applying the global rotation
  /// @param transform The transform to be rotated
  void operator()(Acts::Transform3* transform) const {
    double angleS = randomize != nullptr ? randomize() * angle : angle;
    (*transform) *= Acts::AngleAxis3(angleS, axis);
  }
};

/// This generator applies a local shift
struct LocalShift {
  /// The axis direction
  Acts::AxisDirection axisDirection = Acts::AxisDirection::AxisX;
  /// The shift to be applied
  double shift = 0.0;
  /// Optional randomization function
  std::function<double()> randomize = nullptr;

  /// @brief The call operation applying the local shift
  /// @param transform The transform to be shifted
  void operator()(Acts::Transform3* transform) const {
    double shiftS = randomize != nullptr ? randomize() * shift : shift;

    const auto rotation = transform->rotation();

    switch (axisDirection) {
      case Acts::AxisDirection::AxisX: {
        Acts::Translation3 translation(shiftS * rotation.col(0));
        (*transform) = translation * (*transform);
        break;
      }
      case Acts::AxisDirection::AxisY: {
        Acts::Translation3 translation(shiftS * rotation.col(1));
        (*transform) = translation * (*transform);
        break;
      }
      case Acts::AxisDirection::AxisZ: {
        Acts::Translation3 translation(shiftS * rotation.col(2));
        (*transform) = translation * (*transform);
        break;
      }
      default:
        throw std::runtime_error("LocalShift: Invalid axis direction");
    }
  }
};

}  // namespace ActsExamples::AlignmentGenerator
