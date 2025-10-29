// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/TrackFitting/GsfMixtureReduction.hpp"
#include "Acts/TrackFitting/detail/SymmetricKlDistanceMatrix.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

using namespace Acts;

namespace {

std::tuple<BoundVector, double> computeMeanAndSumOfWeights(
    const std::vector<GsfComponent> &cmps) {
  const auto mean =
      std::accumulate(cmps.begin(), cmps.end(), BoundVector::Zero().eval(),
                      [](auto sum, const auto &cmp) -> BoundVector {
                        return sum + cmp.weight * cmp.boundPars;
                      });

  const double sumOfWeights =
      std::accumulate(cmps.begin(), cmps.end(), 0.0,
                      [](auto sum, const auto &cmp) { return sum + cmp.weight; });

  return std::make_tuple(mean, sumOfWeights);
}

GsfComponent makeDefaultComponent(double weight) {
  GsfComponent cmp;
  cmp.boundPars = BoundVector::Zero();
  cmp.boundCov = BoundSquareMatrix::Identity();
  cmp.weight = weight;
  return cmp;
}

}  // namespace

BOOST_AUTO_TEST_CASE(test_distance_matrix_min_distance) {
  std::vector<GsfComponent> cmps = {
      {1. / 3., BoundVector::Constant(-2.), BoundSquareMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+0.), BoundSquareMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+1.), BoundSquareMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+4.), BoundSquareMatrix::Identity()}};

  const auto proj = [](auto &a) -> decltype(auto) { return a; };
  detail::SymmetricKLDistanceMatrix mat(cmps, proj);

  const auto [i, j] = mat.minDistancePair();
  BOOST_CHECK_EQUAL(std::min(i, j), 1);
  BOOST_CHECK_EQUAL(std::max(i, j), 2);
}

BOOST_AUTO_TEST_CASE(test_distance_matrix_masking) {
  std::vector<GsfComponent> cmps = {
      {1. / 3., BoundVector::Constant(-2.), BoundSquareMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+0.), BoundSquareMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+1.), BoundSquareMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+4.), BoundSquareMatrix::Identity()}};

  const auto proj = [](auto &a) -> decltype(auto) { return a; };
  const std::size_t cmp_to_mask = 2;

  detail::SymmetricKLDistanceMatrix mat_full(cmps, proj);
  mat_full.maskAssociatedDistances(cmp_to_mask);

  cmps.erase(cmps.begin() + cmp_to_mask);
  detail::SymmetricKLDistanceMatrix mat_small(cmps, proj);

  const auto [full_i, full_j] = mat_full.minDistancePair();
  const auto [small_i, small_j] = mat_small.minDistancePair();

  BOOST_CHECK_EQUAL(std::min(full_i, full_j), 0);
  BOOST_CHECK_EQUAL(std::max(full_i, full_j), 1);
  BOOST_CHECK_EQUAL(full_i, small_i);
  BOOST_CHECK_EQUAL(full_j, small_j);
}

BOOST_AUTO_TEST_CASE(test_distance_matrix_recompute_distance) {
  std::vector<GsfComponent> cmps = {
      {1. / 3., BoundVector::Constant(-2.), BoundSquareMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+0.), BoundSquareMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+1.), BoundSquareMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+4.), BoundSquareMatrix::Identity()}};

  const auto proj = [](auto &a) -> decltype(auto) { return a; };
  detail::SymmetricKLDistanceMatrix mat(cmps, proj);

  {
    const auto [i, j] = mat.minDistancePair();
    BOOST_CHECK_EQUAL(std::min(i, j), 1);
    BOOST_CHECK_EQUAL(std::max(i, j), 2);
  }

  cmps[3].boundPars = BoundVector::Constant(0.1);
  mat.recomputeAssociatedDistances(3, cmps, proj);

  {
    const auto [i, j] = mat.minDistancePair();
    BOOST_CHECK_EQUAL(std::min(i, j), 1);
    BOOST_CHECK_EQUAL(std::max(i, j), 3);
  }

  cmps[0].boundPars = BoundVector::Constant(1.01);
  mat.recomputeAssociatedDistances(0, cmps, proj);

  {
    const auto [i, j] = mat.minDistancePair();
    BOOST_CHECK_EQUAL(std::min(i, j), 0);
    BOOST_CHECK_EQUAL(std::max(i, j), 2);
  }
}

BOOST_AUTO_TEST_CASE(test_mixture_reduction) {
  // Assume that the components are on a generic plane surface
  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3{0, 0, 0}, Vector3{1, 0, 0}).planeSurface();
  const std::size_t NComps = 4;
  std::vector<GsfComponent> cmps;

  for (auto i = 0ul; i < NComps; ++i) {
    cmps.push_back(makeDefaultComponent(1.0 / NComps));
  }

  cmps[0].boundPars[eBoundQOverP] = 0.5_GeV;
  cmps[1].boundPars[eBoundQOverP] = 1.5_GeV;
  cmps[2].boundPars[eBoundQOverP] = 3.5_GeV;
  cmps[3].boundPars[eBoundQOverP] = 4.5_GeV;

  // Check start properties
  const auto [mean0, sumOfWeights0] = computeMeanAndSumOfWeights(cmps);

  BOOST_CHECK_CLOSE(mean0[eBoundQOverP], 2.5_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(sumOfWeights0, 1.0, 1.e-8);

  // Reduce by factor of 2 and check if weights and QoP are correct
  reduceMixtureWithKLDistance(cmps, 2, *surface);

  BOOST_CHECK_EQUAL(cmps.size(), 2);

  std::ranges::sort(cmps, {},
                    [](const auto &c) { return c.boundPars[eBoundQOverP]; });
  BOOST_CHECK_CLOSE(cmps[0].boundPars[eBoundQOverP], 1.0_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(cmps[1].boundPars[eBoundQOverP], 4.0_GeV, 1.e-8);

  const auto [mean1, sumOfWeights1] = computeMeanAndSumOfWeights(cmps);

  BOOST_CHECK_CLOSE(mean1[eBoundQOverP], 2.5_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(sumOfWeights1, 1.0, 1.e-8);

  // Reduce by factor of 2 and check if weights and QoP are correct
  reduceMixtureWithKLDistance(cmps, 1, *surface);

  BOOST_CHECK_EQUAL(cmps.size(), 1);
  BOOST_CHECK_CLOSE(cmps[0].boundPars[eBoundQOverP], 2.5_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(cmps[0].weight, 1.0, 1.e-8);
}

BOOST_AUTO_TEST_CASE(test_weight_cut_reduction) {
  std::shared_ptr<PlaneSurface> dummy =
      CurvilinearSurface(Vector3{0, 0, 0}, Vector3{1, 0, 0}).planeSurface();
  std::vector<GsfComponent> cmps;

  // weights do not need to be normalized for this test
  for (auto w : {1.0, 2.0, 3.0, 4.0}) {
    cmps.push_back(makeDefaultComponent(w));
  }

  reduceMixtureLargestWeights(cmps, 2, *dummy);

  BOOST_CHECK_EQUAL(cmps.size(), 2);
  std::ranges::sort(cmps, {}, [](const auto &c) { return c.weight; });

  BOOST_CHECK_EQUAL(cmps[0].weight, 3.0);
  BOOST_CHECK_EQUAL(cmps[1].weight, 4.0);
}

BOOST_AUTO_TEST_CASE(test_naive_vs_optimized) {
  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3{0, 0, 0}, Vector3{1, 0, 0}).planeSurface();
  const std::size_t NComps = 8;

  std::vector<GsfComponent> cmpsOptimized;
  std::vector<GsfComponent> cmpsNaive;

  for (auto i = 0ul; i < NComps; ++i) {
    GsfComponent a = makeDefaultComponent(1.0 / NComps);
    a.boundPars[eBoundQOverP] = static_cast<double>(i) * 1.0_GeV;
    cmpsOptimized.push_back(a);
    cmpsNaive.push_back(a);
  }

  const auto [meanBefore, sumOfWeightsBefore] = computeMeanAndSumOfWeights(cmpsOptimized);

  reduceMixtureWithKLDistance(cmpsOptimized, 2, *surface);
  reduceMixtureWithKLDistanceNaive(cmpsNaive, 2, *surface);

  BOOST_CHECK_EQUAL(cmpsOptimized.size(), 2);
  BOOST_CHECK_EQUAL(cmpsNaive.size(), 2);

  const auto [meanOptimized, sumOfWeightsOptimized] = computeMeanAndSumOfWeights(cmpsOptimized);
  const auto [meanNaive, sumOfWeightsNaive] = computeMeanAndSumOfWeights(cmpsNaive);

  BOOST_CHECK_CLOSE(sumOfWeightsOptimized, sumOfWeightsNaive, 1.e-8);
  BOOST_CHECK_CLOSE(sumOfWeightsOptimized, 1.0, 1.e-8);

  for (auto i = 0ul; i < BoundVector::RowsAtCompileTime; ++i) {
    BOOST_CHECK_CLOSE(meanOptimized[i], meanBefore[i], 1.e-8);
    BOOST_CHECK_CLOSE(meanNaive[i], meanBefore[i], 1.e-8);
  }
}
