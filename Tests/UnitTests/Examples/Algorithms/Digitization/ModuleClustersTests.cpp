// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/BinningData.hpp"
#include "ActsExamples/Digitization/ModuleClusters.hpp"
#include "ActsFatras/Digitization/Channelizer.hpp"

using namespace Acts;
using namespace ActsFatras;
using namespace ActsExamples;

namespace {

auto testDigitizedParameters(bool merge) {
  BinUtility binUtility;
  binUtility +=
      BinningData(BinningOption::open, BinningValue::binX, 20, -10.0f, 10.0f);
  binUtility +=
      BinningData(BinningOption::open, BinningValue::binY, 20, -10.0f, 10.0f);
  std::vector<Acts::BoundIndices> boundIndices = {eBoundLoc0, eBoundLoc1};
  double nsigma = 1;
  bool commonCorner = true;

  ModuleClusters moduleClusters(binUtility, boundIndices, merge, nsigma,
                                commonCorner);

  Cluster cluster;
  cluster.sizeLoc0 = 1;
  cluster.sizeLoc1 = 1;
  cluster.channels = {Cluster::Cell(
      {0, 0}, Channelizer::Segment2D({Vector2(0, 0), Vector2(0, 0)}), 1)};

  DigitizedParameters params;
  params.indices = {eBoundLoc0, eBoundLoc1};
  params.values = {0, 0};
  params.variances = {1, 1};
  params.cluster = {cluster};

  moduleClusters.add(params, 0);
  moduleClusters.add(params, 1);

  return moduleClusters.digitizedParameters();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(DigitizationModuleClustersTests)

BOOST_AUTO_TEST_CASE(digitizedParameters_merging) {
  auto result = testDigitizedParameters(true);
  BOOST_CHECK_EQUAL(result.size(), 1);

  result = testDigitizedParameters(false);
  BOOST_CHECK_EQUAL(result.size(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
