// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/Digitization/ModuleClusters.hpp"
#include "ActsFatras/Digitization/Segmentizer.hpp"

#include <algorithm>
#include <iterator>
#include <memory>

using namespace Acts;
using namespace ActsFatras;
using namespace ActsExamples;

namespace ActsTests {

std::shared_ptr<const IMultiAxis> makeSegmentation() {
  const auto axisX = IAxis::createEquidistant(AxisBoundaryType::Bound, -10.0,
                                              10.0, 20, AxisDirection::AxisX);
  const auto axisY = IAxis::createEquidistant(AxisBoundaryType::Bound, -10.0,
                                              10.0, 20, AxisDirection::AxisY);
  return IMultiAxis::create(*axisX, *axisY);
}

DigitizedParameters makeDigitizationParameters(const Vector2 &position,
                                               const Vector2 &variance,
                                               const IMultiAxis &segmentation) {
  // Cell bins are zero-based while the axis bin indices start at one
  const std::size_t binX = segmentation.getAxis(0).getBin(position.x()) - 1;
  const std::size_t binY = segmentation.getAxis(1).getBin(position.y()) - 1;
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
    double timeVariance, const IMultiAxis &segmentation) {
  DigitizedParameters params =
      makeDigitizationParameters(position, variance, segmentation);
  params.indices.push_back(eBoundTime);
  params.values.push_back(time);
  params.variances.push_back(timeVariance);
  return params;
}

auto testDigitizedParametersWithTwoClusters(bool merge, const Vector2 &firstHit,
                                            const Vector2 &secondHit) {
  std::shared_ptr<const IMultiAxis> segmentation = makeSegmentation();
  std::vector<Acts::BoundIndices> boundIndices = {eBoundLoc0, eBoundLoc1};
  double nsigma = 1;
  bool commonCorner = true;

  ModuleClusters moduleClusters(segmentation, boundIndices, merge, nsigma,
                                commonCorner);

  moduleClusters.add(
      makeDigitizationParameters(firstHit, {1, 1}, *segmentation), 0);
  moduleClusters.add(
      makeDigitizationParameters(secondHit, {1, 1}, *segmentation), 1);

  return moduleClusters.digitizedParameters();
}

auto testDigitizedParametersWithTwoTimedClusters(bool merge,
                                                 const Vector2 &firstHit,
                                                 double firstTime,
                                                 const Vector2 &secondHit,
                                                 double secondTime) {
  std::shared_ptr<const IMultiAxis> segmentation = makeSegmentation();
  std::vector<Acts::BoundIndices> boundIndices = {eBoundLoc0, eBoundLoc1};
  double nsigma = 1;
  bool commonCorner = true;

  ModuleClusters moduleClusters(segmentation, boundIndices, merge, nsigma,
                                commonCorner);

  moduleClusters.add(makeDigitizationParametersWithTime(
                         firstHit, {1, 1}, firstTime, 1, *segmentation),
                     0);
  moduleClusters.add(makeDigitizationParametersWithTime(
                         secondHit, {1, 1}, secondTime, 1, *segmentation),
                     1);

  return moduleClusters.digitizedParameters();
}

BOOST_AUTO_TEST_SUITE(DigitizationSuite)

BOOST_AUTO_TEST_CASE(charge_quantization_response) {
  for (unsigned int bits : {1u, 4u, 8u, 16u, 24u}) {
    const ChargeQuantizationConfig q{bits, 0., 0.25};
    BOOST_CHECK_EQUAL(q.quantize(-1.), 1.);
    BOOST_CHECK_EQUAL(q.quantize(1e9), (1u << bits) - 1u);
  }
  ChargeQuantizationConfig q{4u, 0., 0.25};
  BOOST_CHECK_EQUAL(q.quantize(0.249), 1.);
  BOOST_CHECK_EQUAL(q.quantize(0.25), 2.);
  BOOST_CHECK_EQUAL(q.quantize(0.251), 2.);
  const ChargeQuantizationConfig wide{24u, 0., 0.25};
  for (unsigned int code : {(1u << 24) - 2u, (1u << 24) - 1u}) {
    const double encoded = wide.quantize((code - 1u) * wide.step);
    BOOST_CHECK_EQUAL(encoded, code);
    BOOST_CHECK_EQUAL(static_cast<double>(static_cast<float>(encoded)),
                      encoded);
  }
}

BOOST_AUTO_TEST_CASE(charge_quantization_after_summation_and_before_weighting) {
  const auto segmentation = makeSegmentation();
  const ChargeQuantizationConfig q{4u, 0., 0.25};
  ModuleClusters native(segmentation, {eBoundLoc0, eBoundLoc1}, true, 1., true);
  ModuleClusters encoded(segmentation, {eBoundLoc0, eBoundLoc1}, true, 1., true,
                         q);
  // Two contributions share an address; a third occupies the adjacent channel.
  for (std::size_t i = 0; i < 3; ++i) {
    auto params = makeDigitizationParameters({i < 2 ? 0.5 : 1.5, 0.5}, {1., 1.},
                                             *segmentation);
    params.cluster.channels.front().activation = i < 2 ? 0.3125 : 0.875;
    native.add(params, i);
    encoded.add(params, i);
  }
  const auto reference = native.digitizedParameters();
  const auto result = encoded.digitizedParameters();
  BOOST_REQUIRE_EQUAL(reference.size(), 1);
  BOOST_REQUIRE_EQUAL(result.size(), 1);
  const auto &actual = result.front().first;
  BOOST_REQUIRE_EQUAL(actual.cluster.channels.size(), 2);
  BOOST_REQUIRE_EQUAL(reference.front().first.cluster.channels.size(), 2);
  for (std::size_t i = 0; i < 2; ++i) {
    const auto &cell = reference.front().first.cluster.channels[i];
    BOOST_CHECK(actual.cluster.channels[i].bin == cell.bin);
  }
  // q(0.3125 + 0.3125) = 3, whereas q(0.3125) + q(0.3125) = 4.
  BOOST_CHECK_EQUAL(actual.cluster.channels[0].activation, 3.);
  BOOST_CHECK_EQUAL(actual.cluster.channels[1].activation, 4.);
  BOOST_CHECK_CLOSE(actual.values[0], (3. * 0.5 + 4. * 1.5) / 7., 1e-10);
  BOOST_CHECK_CLOSE(actual.values[1], 0.5, 1e-10);
  BOOST_CHECK_CLOSE(actual.variances[0], 25. / (12. * 49.), 1e-10);
  BOOST_CHECK_CLOSE(actual.variances[1], 25. / (12. * 49.), 1e-10);
  BOOST_CHECK_EQUAL(reference.front().first.cluster.channels[0].activation,
                    0.625);
  BOOST_CHECK_EQUAL(reference.front().first.cluster.channels[1].activation,
                    0.875);
  BOOST_CHECK_CLOSE(reference.front().first.values[0],
                    (0.625 * 0.5 + 0.875 * 1.5) / 1.5, 1e-10);
  BOOST_CHECK(reference.front().second == result.front().second);
}

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
