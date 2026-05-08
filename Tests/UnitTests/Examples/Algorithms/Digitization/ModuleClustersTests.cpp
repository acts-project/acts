// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/Digitization/ModuleClusters.hpp"
#include "ActsFatras/Digitization/Segmentizer.hpp"

#include <algorithm>
#include <iterator>

using namespace Acts;
using namespace ActsFatras;
using namespace ActsExamples;

namespace ActsTests {

DigitizedParameters makeDigitizationParameters(
    const Vector2 &position, const Vector2 &variance,
    const std::vector<DirectedProtoAxis> &segmentation) {
  const std::size_t binX = segmentation[0].bin(position.x());
  const std::size_t binY = segmentation[1].bin(position.y());
  const Segmentizer::Bin2D bin = {
      static_cast<Segmentizer::Bin2D::value_type>(binX),
      static_cast<Segmentizer::Bin2D::value_type>(binY)};
  const Segmentizer::Segment2D segment = {position, position};
  const double activation = 1;
  const Cluster::Cell cell = {bin, segment, activation};

  Cluster cluster;
  cluster.sizeLoc0 = 1;
  cluster.sizeLoc1 = 1;
  cluster.channels = {cell};

  DigitizedParameters params;
  params.indices = {eBoundLoc0, eBoundLoc1};
  params.values = {position.x(), position.y()};
  params.variances = {variance.x(), variance.y()};
  params.cluster = {cluster};

  return params;
}

DigitizedParameters makeDigitizationParametersWithTime(
    const Vector2 &position, const Vector2 &variance, double time,
    double timeVariance, const std::vector<DirectedProtoAxis> &segmentation) {
  DigitizedParameters params =
      makeDigitizationParameters(position, variance, segmentation);
  params.indices.push_back(eBoundTime);
  params.values.push_back(time);
  params.variances.push_back(timeVariance);
  return params;
}

auto testDigitizedParametersWithTwoClusters(bool merge, const Vector2 &firstHit,
                                            const Vector2 &secondHit) {
  std::vector<DirectedProtoAxis> segmentation;
  segmentation.emplace_back(AxisDirection::AxisX, AxisBoundaryType::Open, -10.0,
                            10.0, 20);
  segmentation.emplace_back(AxisDirection::AxisY, AxisBoundaryType::Open, -10.0,
                            10.0, 20);
  std::vector<Acts::BoundIndices> boundIndices = {eBoundLoc0, eBoundLoc1};
  double nsigma = 1;
  bool commonCorner = true;

  ModuleClusters moduleClusters(segmentation, boundIndices, merge, nsigma,
                                commonCorner);

  moduleClusters.add(makeDigitizationParameters(firstHit, {1, 1}, segmentation),
                     0);
  moduleClusters.add(
      makeDigitizationParameters(secondHit, {1, 1}, segmentation), 1);

  return moduleClusters.digitizedParameters();
}

auto testDigitizedParametersWithTwoTimedClusters(bool merge,
                                                 const Vector2 &firstHit,
                                                 double firstTime,
                                                 const Vector2 &secondHit,
                                                 double secondTime) {
  std::vector<DirectedProtoAxis> segmentation;
  segmentation.emplace_back(AxisDirection::AxisX, AxisBoundaryType::Open, -10.0,
                            10.0, 20);
  segmentation.emplace_back(AxisDirection::AxisY, AxisBoundaryType::Open, -10.0,
                            10.0, 20);
  std::vector<Acts::BoundIndices> boundIndices = {eBoundLoc0, eBoundLoc1};
  double nsigma = 1;
  bool commonCorner = true;

  ModuleClusters moduleClusters(segmentation, boundIndices, merge, nsigma,
                                commonCorner);

  moduleClusters.add(makeDigitizationParametersWithTime(
                         firstHit, {1, 1}, firstTime, 1, segmentation),
                     0);
  moduleClusters.add(makeDigitizationParametersWithTime(
                         secondHit, {1, 1}, secondTime, 1, segmentation),
                     1);

  return moduleClusters.digitizedParameters();
}

BOOST_AUTO_TEST_SUITE(DigitizationSuite)

BOOST_AUTO_TEST_CASE(digitizedParameters_merging) {
  // overlapping hits are expected to be merged if turned on
  {
    auto result = testDigitizedParametersWithTwoClusters(true, {0, 0}, {0, 0});
    BOOST_CHECK_EQUAL(result.size(), 1);

    result = testDigitizedParametersWithTwoClusters(false, {0, 0}, {0, 0});
    BOOST_CHECK_EQUAL(result.size(), 2);
  }

  // non overlapping hits are not expected to be merged
  {
    auto result = testDigitizedParametersWithTwoClusters(true, {0, 0}, {5, 0});
    BOOST_CHECK_EQUAL(result.size(), 2);

    result = testDigitizedParametersWithTwoClusters(false, {0, 0}, {5, 0});
    BOOST_CHECK_EQUAL(result.size(), 2);
  }
}

BOOST_AUTO_TEST_CASE(digitizedParameters_merging_with_smeared_time) {
  // hits with geometric (loc0, loc1) and smeared (time) parameters in
  // adjacent cells: compatible times are expected to be merged
  {
    auto result = testDigitizedParametersWithTwoTimedClusters(true, {0, 0}, 1.0,
                                                              {1.05, 0}, 2.0);
    BOOST_REQUIRE_EQUAL(result.size(), 1);

    const auto &[params, simHits] = result.front();
    auto it = std::ranges::find(params.indices, eBoundTime);
    BOOST_REQUIRE(it != params.indices.end());
    auto slot = std::distance(params.indices.begin(), it);
    BOOST_CHECK_CLOSE(params.values.at(slot), 1.5, 1e-6);
    BOOST_CHECK_EQUAL(simHits.size(), 2);
  }

  // incompatible times prevent merging
  {
    auto result = testDigitizedParametersWithTwoTimedClusters(true, {0, 0}, 0.0,
                                                              {1.05, 0}, 10.0);
    BOOST_CHECK_EQUAL(result.size(), 2);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
