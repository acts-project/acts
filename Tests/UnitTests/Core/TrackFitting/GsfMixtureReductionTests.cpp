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
#include "Acts/Utilities/Zip.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {

std::tuple<BoundVector, double> computeMeanAndSumOfWeights(
    const std::vector<GsfComponent> &cmps, const Surface &surface) {
  if (cmps.empty()) {
    return {BoundVector::Zero(), 0.0};
  }

  const double sumOfWeights = std::accumulate(
      cmps.begin(), cmps.end(), 0.0,
      [](auto sum, const auto &cmp) { return sum + cmp.weight; });

  BoundVector mean;
  detail::angleDescriptionSwitch(surface, [&](const auto &desc) {
    auto [meanTmp, cov] = detail::gaussianMixtureMeanCov(
        cmps,
        [](const auto &c) {
          return std::tie(c.weight, c.boundPars, c.boundCov);
        },
        desc);
    mean = meanTmp;
  });

  return std::make_tuple(mean, sumOfWeights);
}

GsfComponent makeDefaultComponent(double weight) {
  GsfComponent cmp;
  cmp.boundPars = BoundVector::Zero();
  cmp.boundCov = BoundSquareMatrix::Identity();
  cmp.weight = weight;
  return cmp;
}

void testReductionEquivalence(const std::vector<GsfComponent> &cmps,
                              std::size_t targetSize, const Surface &surface) {
  std::vector<GsfComponent> cmpsOptimized = cmps;
  std::vector<GsfComponent> cmpsNaive = cmps;

  const auto [meanBefore, sumOfWeightsBefore] =
      computeMeanAndSumOfWeights(cmps, surface);

  reduceMixtureWithKLDistance(cmpsOptimized, targetSize, surface);
  reduceMixtureWithKLDistanceNaive(cmpsNaive, targetSize, surface);

  BOOST_REQUIRE(cmpsOptimized.size() == targetSize);
  BOOST_REQUIRE(cmpsNaive.size() == targetSize);

  // Compare components
  for (const auto &[opt, naiv] : Acts::zip(cmpsOptimized, cmpsNaive)) {
    BOOST_CHECK_CLOSE(opt.weight, naiv.weight, 1.e-8);
    BOOST_CHECK_CLOSE(opt.boundPars[eBoundLoc0], naiv.boundPars[eBoundLoc0],
                      1.e-8);
    BOOST_CHECK_CLOSE(opt.boundPars[eBoundLoc1], naiv.boundPars[eBoundLoc1],
                      1.e-8);
    BOOST_CHECK_CLOSE(opt.boundPars[eBoundPhi], naiv.boundPars[eBoundPhi],
                      1.e-8);
    BOOST_CHECK_CLOSE(opt.boundPars[eBoundTheta], naiv.boundPars[eBoundTheta],
                      1.e-8);
    BOOST_CHECK_CLOSE(opt.boundPars[eBoundQOverP], naiv.boundPars[eBoundQOverP],
                      1.e-8);
    BOOST_CHECK_CLOSE(opt.boundPars[eBoundTime], naiv.boundPars[eBoundTime],
                      1.e-8);
  }

  // Compare mean with mean before
  const auto [meanOptimized, sumOfWeightsOptimized] =
      computeMeanAndSumOfWeights(cmpsOptimized, surface);
  const auto [meanNaive, sumOfWeightsNaive] =
      computeMeanAndSumOfWeights(cmpsNaive, surface);

  BOOST_CHECK_CLOSE(sumOfWeightsOptimized, sumOfWeightsNaive, 1.e-8);
  BOOST_CHECK_CLOSE(sumOfWeightsOptimized, 1.0, 1.e-8);

  BOOST_CHECK_CLOSE(meanOptimized[eBoundLoc0], meanBefore[eBoundLoc0], 1.e-8);
  BOOST_CHECK_CLOSE(meanNaive[eBoundLoc0], meanBefore[eBoundLoc0], 1.e-8);

  BOOST_CHECK_CLOSE(meanOptimized[eBoundLoc1], meanBefore[eBoundLoc1], 1.e-8);
  BOOST_CHECK_CLOSE(meanNaive[eBoundLoc1], meanBefore[eBoundLoc1], 1.e-8);

  // TODO this is concerning, investigate!
  BOOST_CHECK_CLOSE(meanOptimized[eBoundPhi], meanBefore[eBoundPhi], 10);
  BOOST_CHECK_CLOSE(meanNaive[eBoundPhi], meanBefore[eBoundPhi], 10);

  BOOST_CHECK_CLOSE(meanOptimized[eBoundTheta], meanBefore[eBoundTheta], 1.e-8);
  BOOST_CHECK_CLOSE(meanNaive[eBoundTheta], meanBefore[eBoundTheta], 1.e-8);

  BOOST_CHECK_CLOSE(meanOptimized[eBoundQOverP], meanBefore[eBoundQOverP],
                    1.e-8);
  BOOST_CHECK_CLOSE(meanNaive[eBoundQOverP], meanBefore[eBoundQOverP], 1.e-8);

  BOOST_CHECK_CLOSE(meanOptimized[eBoundTime], meanBefore[eBoundTime], 1.e-8);
  BOOST_CHECK_CLOSE(meanNaive[eBoundTime], meanBefore[eBoundTime], 1.e-8);
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

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
  const auto [mean0, sumOfWeights0] =
      computeMeanAndSumOfWeights(cmps, *surface);

  BOOST_CHECK_CLOSE(mean0[eBoundQOverP], 2.5_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(sumOfWeights0, 1.0, 1.e-8);

  // Reduce by factor of 2 and check if weights and QoP are correct
  reduceMixtureWithKLDistance(cmps, 2, *surface);

  BOOST_CHECK_EQUAL(cmps.size(), 2);

  std::ranges::sort(cmps, {},
                    [](const auto &c) { return c.boundPars[eBoundQOverP]; });
  BOOST_CHECK_CLOSE(cmps[0].boundPars[eBoundQOverP], 1.0_GeV, 1.e-8);
  BOOST_CHECK_CLOSE(cmps[1].boundPars[eBoundQOverP], 4.0_GeV, 1.e-8);

  const auto [mean1, sumOfWeights1] =
      computeMeanAndSumOfWeights(cmps, *surface);

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
  const std::size_t NComps = 10;

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> weightDist(0.5, 1.5);
  std::uniform_real_distribution<double> loc0Dist(-10.0, 10.0);
  std::uniform_real_distribution<double> loc1Dist(-10.0, 10.0);
  std::uniform_real_distribution<double> phiDist(-std::numbers::pi,
                                                 std::numbers::pi);
  std::uniform_real_distribution<double> thetaDist(0.0, std::numbers::pi);
  std::uniform_real_distribution<double> qopDist(0.1, 5.0);

  std::vector<GsfComponent> cmps;
  double weightSum = 0.0;

  for (auto i = 0ul; i < NComps; ++i) {
    GsfComponent cmp = makeDefaultComponent(weightDist(rng));
    cmp.boundPars[eBoundLoc0] = loc0Dist(rng);
    cmp.boundPars[eBoundLoc1] = loc1Dist(rng);
    cmp.boundPars[eBoundPhi] = phiDist(rng);
    cmp.boundPars[eBoundTheta] = thetaDist(rng);
    cmp.boundPars[eBoundQOverP] = qopDist(rng);
    cmp.boundPars[eBoundTime] = 0.0;
    weightSum += cmp.weight;
    cmps.push_back(cmp);
  }

  for (auto &cmp : cmps) {
    cmp.weight /= weightSum;
  }

  for (std::size_t targetSize = 9; targetSize > 0; --targetSize) {
    testReductionEquivalence(cmps, targetSize, *surface);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
