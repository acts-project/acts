// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
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

BOOST_AUTO_TEST_SUITE(DigitizationModuleClustersTests)

BOOST_AUTO_TEST_CASE(merge) {
  BinUtility binUtility;
  binUtility += BinningData(binX, -10.0f, 10.0f);
  binUtility += BinningData(binY, -10.0f, 10.0f);
  std::vector<Acts::BoundIndices> boundIndices = {eBoundLoc0, eBoundLoc1};
  bool merge = true;
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

  auto result = moduleClusters.digitizedParameters();

  BOOST_CHECK_EQUAL(result.size(), 1);
}

BOOST_AUTO_TEST_CASE(nomerge) {
  BinUtility binUtility;
  binUtility += BinningData(binX, -10.0f, 10.0f);
  binUtility += BinningData(binY, -10.0f, 10.0f);
  std::vector<Acts::BoundIndices> boundIndices = {eBoundLoc0, eBoundLoc1};
  bool merge = false;
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

  auto result = moduleClusters.digitizedParameters();

  BOOST_CHECK_EQUAL(result.size(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
