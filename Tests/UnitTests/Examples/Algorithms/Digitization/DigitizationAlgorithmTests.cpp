// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"

using namespace Acts;
using namespace ActsExamples;

BOOST_AUTO_TEST_SUITE(DigitizationAlgorithmTests)

BOOST_AUTO_TEST_CASE(merge) {
  // constants
  const auto simhitHandle = "simhits";
  const auto clustersHandle = "clusters";

  // setup sequencer chain
  bool merge = true;
  double sigma = 1;
  bool commonCorner = true;
  DigitizationConfig config(merge, sigma, commonCorner);
  DigitizationAlgorithm algorithm(config, Logging::Level::VERBOSE);

  // sequencer init run
  algorithm.initialize();

  // sequencer first event
  WhiteBoard store;
  AlgorithmContext ctx(0, 0, store);
  // sequencer runs our reader
  {
    WriteDataHandle<SimHitContainer> simHitContainerWriteHandle(
        nullptr, "SimHitContainer");
    simHitContainerWriteHandle.initialize(simhitHandle);

    SimHitContainer simHits;
    simHits.emplace(Acts::GeometryIdentifier(), SimBarcode(), Vector4::Zero(),
                    Vector4::Zero(), Vector4::Zero());
    simHits.emplace(Acts::GeometryIdentifier(), SimBarcode(), Vector4::Zero(),
                    Vector4::Zero(), Vector4::Zero());
    simHitContainerWriteHandle(ctx, std::move(simHits));
  }
  // sequencer runs the geometric digitization
  algorithm.execute(ctx);
  // sequencer runs our writer
  {
    ReadDataHandle<SimHitContainer> clusterReadHandle(nullptr, "Clusters");
    clusterReadHandle.initialize(clustersHandle);

    auto clusters = clusterReadHandle(ctx);

    BOOST_CHECK_EQUAL(clusters.size(), 1);
  }
}

BOOST_AUTO_TEST_SUITE_END()
