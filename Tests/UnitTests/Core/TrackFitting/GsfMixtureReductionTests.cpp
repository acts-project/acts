// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/TrackFitting/detail/KLMixtureReduction.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

struct DummyComponent {
  double weight = 0.0;
  BoundVector boundPars = BoundVector::Zero();
  BoundSymMatrix boundCov = BoundSymMatrix::Identity();
};

BOOST_AUTO_TEST_CASE(test_distance_matrix_min_distance) {
  std::vector<DummyComponent> cmps = {
      {1. / 3., BoundVector::Constant(-2.), BoundSymMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+0.), BoundSymMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+1.), BoundSymMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+4.), BoundSymMatrix::Identity()}};

  const auto proj = [](auto &a) -> decltype(auto) { return a; };
  detail::SymmetricKLDistanceMatrix mat(cmps, proj);

  const auto [i, j] = mat.minDistancePair();
  BOOST_CHECK(std::min(i, j) == 1);
  BOOST_CHECK(std::max(i, j) == 2);
}

BOOST_AUTO_TEST_CASE(test_distance_matrix_masking) {
  std::vector<DummyComponent> cmps = {
      {1. / 3., BoundVector::Constant(-2.), BoundSymMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+0.), BoundSymMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+1.), BoundSymMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+4.), BoundSymMatrix::Identity()}};

  const auto proj = [](auto &a) -> decltype(auto) { return a; };
  const std::size_t cmp_to_mask = 2;

  detail::SymmetricKLDistanceMatrix mat_full(cmps, proj);
  mat_full.maskAssociatedDistances(cmp_to_mask);

  cmps.erase(cmps.begin() + cmp_to_mask);
  detail::SymmetricKLDistanceMatrix mat_small(cmps, proj);

  const auto [full_i, full_j] = mat_full.minDistancePair();
  const auto [small_i, small_j] = mat_small.minDistancePair();

  BOOST_CHECK(std::min(full_i, full_j) == 0);
  BOOST_CHECK(std::max(full_i, full_j) == 1);
  BOOST_CHECK(full_i == small_i);
  BOOST_CHECK(full_j == small_j);
}

BOOST_AUTO_TEST_CASE(test_distance_matrix_recompute_distance) {
  std::vector<DummyComponent> cmps = {
      {1. / 3., BoundVector::Constant(-2.), BoundSymMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+0.), BoundSymMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+1.), BoundSymMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+4.), BoundSymMatrix::Identity()}};

  const auto proj = [](auto &a) -> decltype(auto) { return a; };
  detail::SymmetricKLDistanceMatrix mat(cmps, proj);

  {
    const auto [i, j] = mat.minDistancePair();
    BOOST_CHECK(std::min(i, j) == 1);
    BOOST_CHECK(std::max(i, j) == 2);
  }

  cmps[3].boundPars = BoundVector::Constant(0.1);
  mat.recomputeAssociatedDistances(3, cmps, proj);

  {
    const auto [i, j] = mat.minDistancePair();
    BOOST_CHECK(std::min(i, j) == 1);
    BOOST_CHECK(std::max(i, j) == 3);
  }

  cmps[0].boundPars = BoundVector::Constant(1.01);
  mat.recomputeAssociatedDistances(0, cmps, proj);

  {
    const auto [i, j] = mat.minDistancePair();
    BOOST_CHECK(std::min(i, j) == 0);
    BOOST_CHECK(std::max(i, j) == 2);
  }
}

BOOST_AUTO_TEST_CASE(test_mixture_reduction) {
  auto meanAndSumOfWeights = [](const auto &cmps) {
    const auto mean = std::accumulate(
        cmps.begin(), cmps.end(), Acts::BoundVector::Zero().eval(),
        [](auto sum, const auto &cmp) -> Acts::BoundVector {
          return sum + cmp.weight * cmp.boundPars;
        });

    const double sumOfWeights = std::accumulate(
        cmps.begin(), cmps.end(), 0.0,
        [](auto sum, const auto &cmp) { return sum + cmp.weight; });

    return std::make_tuple(mean, sumOfWeights);
  };

  // Do not bother with circular angles in this test
  const auto desc = std::tuple<>{};

  // Need this projection, since we need to write to the lvalue references which
  // isn't possible through Identity / std::identity due to perfect forwarding
  const auto proj = [](auto &a) -> decltype(auto) { return a; };

  const std::size_t NComps = 4;
  std::vector<DummyComponent> cmps;

  for (auto i = 0ul; i < NComps; ++i) {
    DummyComponent a;
    a.boundPars = BoundVector::Zero();
    a.boundCov = BoundSymMatrix::Identity();
    a.weight = 1.0 / NComps;
    cmps.push_back(a);
  }

  cmps[0].boundPars[eBoundQOverP] = 0.5_GeV;
  cmps[1].boundPars[eBoundQOverP] = 1.5_GeV;
  cmps[2].boundPars[eBoundQOverP] = 3.5_GeV;
  cmps[3].boundPars[eBoundQOverP] = 4.5_GeV;

  // Check start properties
  const auto [mean0, sumOfWeights0] = meanAndSumOfWeights(cmps);

  BOOST_CHECK_CLOSE(mean0[eBoundQOverP], 2.5_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(sumOfWeights0, 1.0, 1.e-8);

  // Reduce by factor of 2 and check if weights and QoP are correct
  Acts::detail::reduceWithKLDistance(cmps, 2, proj, desc);

  BOOST_CHECK(cmps.size() == 2);

  std::sort(cmps.begin(), cmps.end(), [](const auto &a, const auto &b) {
    return a.boundPars[eBoundQOverP] < b.boundPars[eBoundQOverP];
  });
  BOOST_CHECK_CLOSE(cmps[0].boundPars[eBoundQOverP], 1.0_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(cmps[1].boundPars[eBoundQOverP], 4.0_GeV, 1.e-8);

  const auto [mean1, sumOfWeights1] = meanAndSumOfWeights(cmps);

  BOOST_CHECK_CLOSE(mean1[eBoundQOverP], 2.5_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(sumOfWeights1, 1.0, 1.e-8);

  // Reduce by factor of 2 and check if weights and QoP are correct
  Acts::detail::reduceWithKLDistance(cmps, 1, proj, desc);

  BOOST_CHECK(cmps.size() == 1);
  BOOST_CHECK_CLOSE(cmps[0].boundPars[eBoundQOverP], 2.5_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(cmps[0].weight, 1.0, 1.e-8);
}
