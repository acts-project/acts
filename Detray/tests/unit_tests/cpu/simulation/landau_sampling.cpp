// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/test/utils/landau_distribution.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <iterator>

using namespace detray;

/// Test the landau distribution
/// Validate that the mpv = -0.22278 for mu = 0 and sigma=1

// Test class for covariance transport
template <typename T>
class detray_simulation_LandauSamplingValidation : public ::testing::Test {
 public:
  using scalar_type = T;

  // Function for getting the corresponding index of a value
  std::size_t get_index(const scalar_type value) const {
    return static_cast<std::size_t>((value - min) / bin_size);
  }

  // Landau distribution with (mu = 0, sigma = 1) has the most probable value
  // of -0.22278
  constexpr static const scalar_type mu = 0.f;
  constexpr static const scalar_type sigma = 1.f;
  constexpr static const scalar_type mpv = -0.22278f;

  // Binning information for counting
  constexpr static const double bin_size = 0.05;
  constexpr static const double min = -2.;
  constexpr static const double max = 2.;
  constexpr static const std::size_t n_bins =
      std::size_t((max - min) / bin_size);
};

// Test for float and double types
using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(detray_simulation_LandauSamplingValidation, TestTypes, );

TYPED_TEST(detray_simulation_LandauSamplingValidation, landau_sampling) {
  // Random generator
  std::random_device rd{};
  std::mt19937_64 generator{rd()};
  generator.seed(0u);

  // Landau distribution for sampling
  landau_distribution<typename TestFixture::scalar_type> ld;

  // Make sure that the number of n_bins is 2 - (-2) / 0.05 = 80
  EXPECT_EQ(this->n_bins, 80u);

  // Counter vector
  std::vector<int> counter(this->n_bins, 0);

  // Sampling and counting
  std::size_t n_samples = 10000000u;
  const auto minf = static_cast<typename TestFixture::scalar_type>(this->min);
  const auto maxf = static_cast<typename TestFixture::scalar_type>(this->max);
  for (std::size_t i = 0u; i < n_samples; i++) {
    const auto sa = ld(generator, this->mu, this->sigma);

    if (sa > minf && sa < maxf) {
      const std::size_t index = this->get_index(sa);
      counter[index]++;
    }
  }

  const std::size_t mpv_index = this->get_index(this->mpv);

  const auto max_index = static_cast<std::size_t>(
      std::distance(counter.begin(), std::ranges::max_element(counter)));

  // Bin range for i index: [ -2 + 0.05 * i, -2 + 0.05 * (i+1) ]
  // Bin range for i = 35 : [ -0.25, -0.2] which includes mpv (-0.22278)
  EXPECT_EQ(mpv_index, 35u);
  // Make sure that max and mpv index is close to each other
  EXPECT_TRUE(max_index == mpv_index || max_index == mpv_index - 1u);
}
