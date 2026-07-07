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
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/TrackFitting/GsfMixtureReduction.hpp"
#include "Acts/TrackFitting/detail/GsfComponentMerging.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <numbers>
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
  detail::Gsf::angleDescriptionSwitch(surface, [&](const auto &desc) {
    auto [meanTmp, cov] = detail::Gsf::mergeGaussianMixtureMeanCov(
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
  cmp.boundCov = BoundMatrix::Identity();
  cmp.weight = weight;
  return cmp;
}

void testReductionEquivalence(std::vector<GsfComponent> &cmpsRef,
                              std::vector<GsfComponent> &cmpsTest) {
  BOOST_REQUIRE(cmpsRef.size() == cmpsTest.size());

  // sort by weight since the order of components is not guaranteed to be the
  // same after reduction
  std::ranges::sort(cmpsRef, {}, [](const auto &c) { return c.weight; });
  std::ranges::sort(cmpsTest, {}, [](const auto &c) { return c.weight; });

  constexpr static double tol = 1.e-8;

  // Compare components
  for (const auto &[ref, test] : Acts::zip(cmpsRef, cmpsTest)) {
    BOOST_CHECK_CLOSE(ref.weight, test.weight, tol);
    BOOST_CHECK_CLOSE(ref.boundPars[eBoundLoc0], test.boundPars[eBoundLoc0],
                      tol);
    BOOST_CHECK_CLOSE(ref.boundPars[eBoundLoc1], test.boundPars[eBoundLoc1],
                      tol);
    BOOST_CHECK_CLOSE(ref.boundPars[eBoundPhi], test.boundPars[eBoundPhi], tol);
    BOOST_CHECK_CLOSE(ref.boundPars[eBoundTheta], test.boundPars[eBoundTheta],
                      tol);
    BOOST_CHECK_CLOSE(ref.boundPars[eBoundQOverP], test.boundPars[eBoundQOverP],
                      tol);
    BOOST_CHECK_CLOSE(ref.boundPars[eBoundTime], test.boundPars[eBoundTime],
                      tol);
  }
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

BOOST_AUTO_TEST_CASE(test_distance_matrix_min_distance) {
  std::vector<GsfComponent> cmps = {
      {1. / 3., BoundVector::Constant(-2.), BoundMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+0.), BoundMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+1.), BoundMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+4.), BoundMatrix::Identity()}};

  detail::Gsf::SymmetricKLDistanceMatrix mat(cmps);

  const auto [i, j] = mat.minDistancePair();
  BOOST_CHECK_EQUAL(std::min(i, j), 1);
  BOOST_CHECK_EQUAL(std::max(i, j), 2);
}

BOOST_AUTO_TEST_CASE(test_distance_matrix_masking) {
  std::vector<GsfComponent> cmps = {
      {1. / 3., BoundVector::Constant(-2.), BoundMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+0.), BoundMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+1.), BoundMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+4.), BoundMatrix::Identity()}};

  const std::size_t cmp_to_mask = 2;

  detail::Gsf::SymmetricKLDistanceMatrix mat_full(cmps);
  mat_full.maskAssociatedDistances(cmp_to_mask);

  cmps.erase(cmps.begin() + cmp_to_mask);
  detail::Gsf::SymmetricKLDistanceMatrix mat_small(cmps);

  const auto [full_i, full_j] = mat_full.minDistancePair();
  const auto [small_i, small_j] = mat_small.minDistancePair();

  BOOST_CHECK_EQUAL(std::min(full_i, full_j), 0);
  BOOST_CHECK_EQUAL(std::max(full_i, full_j), 1);
  BOOST_CHECK_EQUAL(full_i, small_i);
  BOOST_CHECK_EQUAL(full_j, small_j);
}

BOOST_AUTO_TEST_CASE(test_distance_matrix_recompute_distance) {
  std::vector<GsfComponent> cmps = {
      {1. / 3., BoundVector::Constant(-2.), BoundMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+0.), BoundMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+1.), BoundMatrix::Identity()},
      {1. / 3., BoundVector::Constant(+4.), BoundMatrix::Identity()}};

  detail::Gsf::SymmetricKLDistanceMatrix mat(cmps);

  {
    const auto [i, j] = mat.minDistancePair();
    BOOST_CHECK_EQUAL(std::min(i, j), 1);
    BOOST_CHECK_EQUAL(std::max(i, j), 2);
  }

  cmps[3].boundPars = BoundVector::Constant(0.1);
  mat.recomputeAssociatedDistances(3, cmps);

  {
    const auto [i, j] = mat.minDistancePair();
    BOOST_CHECK_EQUAL(std::min(i, j), 1);
    BOOST_CHECK_EQUAL(std::max(i, j), 3);
  }

  cmps[0].boundPars = BoundVector::Constant(1.01);
  mat.recomputeAssociatedDistances(0, cmps);

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
    GsfComponent a;
    a.boundPars = BoundVector::Zero();
    a.boundCov = BoundMatrix::Identity();
    a.weight = 1.0 / NComps;
    cmps.push_back(a);
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
    GsfComponent a;
    a.boundPars = BoundVector::Zero();
    a.boundCov = BoundMatrix::Identity();
    a.weight = w;
    cmps.push_back(a);
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
    std::vector<GsfComponent> cmpsOptimized = cmps;
    std::vector<GsfComponent> cmpsNaive = cmps;

    reduceMixtureWithKLDistance(cmpsOptimized, targetSize, *surface);
    reduceMixtureWithKLDistanceNaive(cmpsNaive, targetSize, *surface);
    BOOST_CHECK_EQUAL(cmpsNaive.size(), targetSize);

    testReductionEquivalence(cmpsNaive, cmpsOptimized);
  }
}

BOOST_AUTO_TEST_CASE(test_exact_tie_breaking) {
  // Four exactly identical components: every pairwise KL distance is exactly
  // 0, so every merge step is a tie. The optimized reducer does not attempt
  // to match naive's specific tie-break convention (which pair is merged
  // first among exact ties is otherwise unspecified) -- but reducing to a
  // single component here still must equal naive's result regardless of
  // merge order, since merging any subset of exactly-identical Gaussians is
  // a commutative, associative operation (the weighted mean/covariance of
  // identical components doesn't depend on pairing order). This is a
  // regression guard for the underlying merge math, not for tie-break
  // ordering.
  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3{0, 0, 0}, Vector3{1, 0, 0}).planeSurface();

  std::vector<GsfComponent> cmps;
  for (int i = 0; i < 4; ++i) {
    GsfComponent cmp = makeDefaultComponent(0.25);
    cmp.boundPars[eBoundQOverP] = 1.0;
    cmp.boundCov(eBoundQOverP, eBoundQOverP) = 1.0;
    cmps.push_back(cmp);
  }

  std::vector<GsfComponent> cmpsNaive = cmps;
  reduceMixtureWithKLDistanceNaive(cmpsNaive, 1, *surface);

  std::vector<GsfComponent> cmpsOptimized = cmps;
  reduceMixtureWithKLDistance(cmpsOptimized, 1, *surface);

  testReductionEquivalence(cmpsNaive, cmpsOptimized);
}

BOOST_AUTO_TEST_CASE(test_naive_vs_optimized_stress) {
  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3{0, 0, 0}, Vector3{1, 0, 0}).planeSurface();

  std::mt19937 rng(1234);
  std::uniform_real_distribution<double> weightDist(0.5, 1.5);
  std::uniform_real_distribution<double> qopDist(0.1, 5.0);
  std::uniform_real_distribution<double> covDist(0.5, 2.0);

  for (int seed = 0; seed < 500; ++seed) {
    const std::size_t NComps = 72;
    std::vector<GsfComponent> cmps;
    double weightSum = 0.0;

    for (auto i = 0ul; i < NComps; ++i) {
      GsfComponent cmp = makeDefaultComponent(weightDist(rng));
      // Force frequent exact ties/duplicates by quantizing qop and cov, and
      // by making every third component an exact duplicate of a previous
      // one, similar to how real GSF material-mixture components can be
      // near/exactly-degenerate.
      double qop = qopDist(rng);
      if (i % 3 == 0 && !cmps.empty()) {
        cmp.boundPars[eBoundQOverP] = cmps[i / 2].boundPars[eBoundQOverP];
        cmp.boundCov(eBoundQOverP, eBoundQOverP) =
            cmps[i / 2].boundCov(eBoundQOverP, eBoundQOverP);
      } else {
        cmp.boundPars[eBoundQOverP] = qop;
        cmp.boundCov(eBoundQOverP, eBoundQOverP) = covDist(rng);
      }
      weightSum += cmp.weight;
      cmps.push_back(cmp);
    }

    for (auto &cmp : cmps) {
      cmp.weight /= weightSum;
    }

    for (std::size_t targetSize : {1ul, 2ul, 12ul, 40ul, 71ul}) {
      std::vector<GsfComponent> cmpsOptimized = cmps;
      std::vector<GsfComponent> cmpsNaive = cmps;

      reduceMixtureWithKLDistance(cmpsOptimized, targetSize, *surface);
      reduceMixtureWithKLDistanceNaive(cmpsNaive, targetSize, *surface);
      BOOST_REQUIRE_EQUAL(cmpsNaive.size(), targetSize);

      testReductionEquivalence(cmpsNaive, cmpsOptimized);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_real_odd_mixture_repro) {
  // Exact 37-component mixture captured from a real GSF fit on the
  // OpenDataDetector where a naive-vs-optimized crosscheck initially flagged
  // a "mismatch". A tie-resistant comparison (sorting by (weight, qop)
  // rather than weight alone) shows the two results are actually identical
  // as sets -- two surviving components share the exact same weight but
  // differ (by 1 ULP) in q/p, so a weight-only sort doesn't produce a
  // canonical order and a naive row-by-row comparison sees a false
  // "mismatch" that is actually just those two rows swapped.
  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3{0, 0, 0}, Vector3{1, 0, 0}).planeSurface();

  const std::vector<std::array<double, 3>> raw = {
      {0.00012655666579041247, 0.25417488769818658, 0.0020918704314837823},
      {0.00015968851779559243, 0.30629702733095354, 0.0020538560685349071},
      {0.00029661316560717379, 0.12812293547922468, 0.00056182609781417687},
      {0.00033688100042926315, 0.11084709577220192, 0.00029893983583882753},
      {0.0010318941615459665, 0.028859063193333447, 7.1302719651786278e-07},
      {0.00014254500574071701, 0.40901951916871326, 0.0032869075295659045},
      {0.00016189674800663739, 0.35386815387761877, 0.00060301272524523138},
      {0.00025420415620657226, 0.33852496378268576, 0.00039097333074605812},
      {0.002927468530316577, 0.33552622385425629, 0.00038042761330121718},
      {0.0001375461242792763, 0.2543563861886427, 0.0020978048143056312},
      {0.00032236935988056527, 0.12821442413726197, 0.00056288151395111809},
      {0.00036613382363523585, 0.11092624828304315, 0.00029915658855550465},
      {0.00057488940087665997, 0.10611665224849191, 0.00027832117183572614},
      {0.0066205342007456735, 0.10517664381500796, 0.0002772849284088526},
      {0.00028884025559003482, 0.17780801717105207, 0.0011000045525178595},
      {0.00067696031842999201, 0.089628387621470684, 0.00034993058928141927},
      {0.00076886365477235098, 0.077543076983460302, 0.00022105561567908419},
      {0.0012072406580476455, 0.074180925974204129, 0.00021087393071708628},
      {0.013902837340240878, 0.073523812530749133, 0.00021036754758495547},
      {0.00014254500574071701, 0.4090195191687131, 0.00056154005300243894},
      {0.00067696031842999201, 0.08962838762147067, 0.00030839088704642868},
      {0.0015866045811124363, 0.045179334404795594, 0.00011780419161756486},
      {0.0018020001523011235, 0.039087444266076778, 8.5058280485853386e-05},
      {0.002829432547842417, 0.037392671562944742, 8.2471210695519588e-05},
      {0.032584340343086715, 0.037061437801057591, 8.2342543534007535e-05},
      {0.00016189674800663739, 0.35386815387761866, 0.00041902767216965071},
      {0.00076886365477235087, 0.077543076983460302, 0.00022954417399631549},
      {0.0018020001523011235, 0.039087444266076785, 8.6889018018308146e-05},
      {0.0020466375728075349, 0.033816972281280139, 6.2378529237207205e-05},
      {0.0032135529815267004, 0.032350719304115642, 6.0442093728760039e-05},
      {0.037007952050475558, 0.032064148433221348, 6.0345785673502317e-05},
      {0.0031816726865231495, 0.33576581234380087, 0.00038145716959681844},
      {0.015110077998288523, 0.073576313521543396, 0.00021086410086757797},
      {0.035413772890929125, 0.037087902182294651, 8.2430854368852665e-05},
      {0.04022150503200226, 0.032087044410779829, 6.036392222156534e-05},
      {0.063154287370857834, 0.030695798500164729, 5.8620538343988411e-05},
      {0.72729805677337689, 0.030423887337191948, 5.853383165528916e-05},
  };

  std::vector<GsfComponent> cmps;
  for (const auto &[weight, qop, cov] : raw) {
    GsfComponent cmp = makeDefaultComponent(weight);
    cmp.boundPars[eBoundQOverP] = qop;
    cmp.boundCov(eBoundQOverP, eBoundQOverP) = cov;
    cmps.push_back(cmp);
  }

  std::vector<GsfComponent> cmpsNaive = cmps;
  std::vector<GsfComponent> cmpsOptimized = cmps;
  reduceMixtureWithKLDistanceNaive(cmpsNaive, 12, *surface);
  reduceMixtureWithKLDistance(cmpsOptimized, 12, *surface);

  BOOST_REQUIRE_EQUAL(cmpsNaive.size(), cmpsOptimized.size());

  const auto sortKey = [](const auto &c) {
    return std::make_pair(c.weight, c.boundPars[eBoundQOverP]);
  };
  std::ranges::sort(cmpsNaive, {}, sortKey);
  std::ranges::sort(cmpsOptimized, {}, sortKey);

  constexpr double tol = 1.e-8;
  for (const auto &[ref, test] : Acts::zip(cmpsNaive, cmpsOptimized)) {
    BOOST_CHECK_CLOSE(ref.weight, test.weight, tol);
    BOOST_CHECK_CLOSE(ref.boundPars[eBoundQOverP], test.boundPars[eBoundQOverP],
                      tol);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
