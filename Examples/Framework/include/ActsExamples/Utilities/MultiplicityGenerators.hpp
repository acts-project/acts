// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <random>

namespace ActsExamples {

/// @brief Generator interface for event multiplicity of vertices
struct MultiplicityGenerator {
  /// @brief Virtual destructor required
  virtual ~MultiplicityGenerator() = default;
  /// @brief Generate the multiplicity of vertices
  ///
  /// @param rng Shared random number generator instance
  /// @return std::size_t The multiplicity for the event
  virtual std::size_t operator()(RandomEngine& rng) const = 0;
};

struct FixedMultiplicityGenerator : public MultiplicityGenerator {
  std::size_t n = 1;

  explicit FixedMultiplicityGenerator(std::size_t _n) : n{_n} {}
  FixedMultiplicityGenerator() = default;

  std::size_t operator()(RandomEngine& /*rng*/) const override { return n; }
};

struct PoissonMultiplicityGenerator : public MultiplicityGenerator {
  double mean = 1;
  explicit PoissonMultiplicityGenerator(double _mean) : mean{_mean} {}
  PoissonMultiplicityGenerator() = default;

  std::size_t operator()(RandomEngine& rng) const override {
    return (0 < mean) ? std::poisson_distribution<std::size_t>(mean)(rng) : 0;
  }
};

}  // namespace ActsExamples
