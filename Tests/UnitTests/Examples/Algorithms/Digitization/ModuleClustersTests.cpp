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
#include "ActsExamples/Digitization/ModuleClusters.hpp"
#include "ActsFatras/Digitization/Segmentizer.hpp"

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

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
