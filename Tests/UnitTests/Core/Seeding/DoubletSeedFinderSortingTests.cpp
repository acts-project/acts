// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/DoubletSeedFinder.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <random>

BOOST_AUTO_TEST_CASE(IntegerDoubletSortPreservesReferenceOrder) {
  std::vector<std::uint64_t> scratch;
  std::default_random_engine rng(42571);
  std::uniform_real_distribution<float> dist(-30.f, 30.f);
  auto check = [&](const std::vector<float>& values, bool interior = false) {
    Acts::DoubletsForMiddleSp doublets;
    for (std::size_t i = 0; i < values.size(); ++i) {
      doublets.emplace_back(i, values[i], 1.f, 0.f, 0.f, 0.f, 0.f, 0.f);
    }
    const Acts::DoubletsForMiddleSp::IndexRange range{
        interior ? 1u : 0u,
        static_cast<unsigned int>(values.size() - (interior ? 1u : 0u))};
    std::vector<Acts::DoubletsForMiddleSp::IndexAndCotTheta> reference,
        candidate;
    doublets.sortByCotTheta(range, reference);
    doublets.sortByCotTheta(range, candidate, scratch);
    BOOST_REQUIRE_EQUAL(reference.size(), candidate.size());
    for (std::size_t i = 0; i < reference.size(); ++i) {
      BOOST_CHECK_EQUAL(reference[i].index, candidate[i].index);
      BOOST_CHECK_EQUAL(std::bit_cast<std::uint32_t>(reference[i].cotTheta),
                        std::bit_cast<std::uint32_t>(candidate[i].cotTheta));
    }
  };
  for (unsigned int size : {0u, 1u, 2u, 31u, 32u, 128u, 2048u}) {
    std::vector<float> values(size);
    std::ranges::generate(values, [&] { return dist(rng); });
    check(values);
    if (size > 2u) {
      check(values, true);
    }
    std::ranges::sort(values);
    check(values);
    std::ranges::reverse(values);
    check(values);
  }
  check({1.f, 1.f, 1.f, -1.f, -1.f, 1.f});
  check({-0.f, 0.f, -1.f, 1.f, -0.f, 0.f});
  check({std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
         std::numeric_limits<float>::denorm_min(),
         -std::numeric_limits<float>::denorm_min()});
  check({std::numeric_limits<float>::infinity(), 1.f,
         -std::numeric_limits<float>::infinity(), 0.f});
}
