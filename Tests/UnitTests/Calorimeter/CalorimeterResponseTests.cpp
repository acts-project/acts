// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "ActsCalorimeter/Digitization/CalorimeterResponse.hpp"
#include "ActsExamples/Calorimeter/CalorimeterDigitizationAlgorithm.hpp"
#include "ActsTests/CommonHelpers/WhiteBoardUtilities.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ActsTests {
namespace {
using namespace Acts::UnitLiterals;
using namespace ActsCalorimeter;

SimCalorimeterHit deposit(std::uint64_t cellId, double energy, double time) {
  return {cellId, Acts::Vector3{1500_mm, 0, 0}, energy, time};
}
}  // namespace

BOOST_AUTO_TEST_SUITE(CalorimeterResponseSuite)

BOOST_AUTO_TEST_CASE(MergeBeforeThresholdAndCalibrateOnce) {
  CalorimeterResponse::Config config;
  config.energyScale = 2;
  config.energyThreshold = 0.05_GeV;
  const CalorimeterResponse response(config);
  const SimCalorimeterHitContainer input{deposit(42, 0.01_GeV, 1_ns),
                                         deposit(42, 0.02_GeV, 4_ns),
                                         deposit(43, 0.001_GeV, 2_ns)};
  const auto output = response(input);
  BOOST_REQUIRE_EQUAL(output.size(), 1);
  BOOST_CHECK_EQUAL(output[0].cellId, 42);
  BOOST_CHECK_CLOSE(output[0].energy / 1_GeV, 0.06, 1e-10);
  BOOST_CHECK_CLOSE(output[0].time / 1_ns, 3, 1e-10);
  BOOST_CHECK(output[0].position == input[0].position);
  const std::vector<std::size_t> sources{0, 1};
  BOOST_CHECK_EQUAL_COLLECTIONS(output[0].sourceIndices.begin(),
                                output[0].sourceIndices.end(), sources.begin(),
                                sources.end());
  BOOST_CHECK_EQUAL(input[0].depositedEnergy, 0.01_GeV);
}

BOOST_AUTO_TEST_CASE(TimeWindowAndThresholdBoundaries) {
  CalorimeterResponse::Config config;
  config.timeMin = 0_ns;
  config.timeMax = 10_ns;
  config.energyThreshold = 0.02_GeV;
  const CalorimeterResponse response(config);
  const SimCalorimeterHitContainer input{
      deposit(7, 0.01_GeV, 0_ns), deposit(7, 0.01_GeV, 10_ns),
      deposit(7, 100_GeV, -1_ns), deposit(7, 100_GeV, 11_ns),
      deposit(8, 0_GeV, 5_ns)};
  const auto output = response(input);
  BOOST_REQUIRE_EQUAL(output.size(), 1);
  BOOST_CHECK_EQUAL(output[0].energy, config.energyThreshold);
  BOOST_CHECK_CLOSE(output[0].time / 1_ns, 5, 1e-10);
  BOOST_CHECK_EQUAL(output[0].sourceIndices.size(), 2);
  BOOST_CHECK(response(SimCalorimeterHitContainer{}).empty());
}

BOOST_AUTO_TEST_CASE(FullCellIdsAndEventIsolation) {
  const auto highId = std::numeric_limits<std::uint64_t>::max();
  const CalorimeterResponse response(CalorimeterResponse::Config{});
  SimCalorimeterHitContainer input{deposit(highId, 2_GeV, 0_ns),
                                   deposit(0, 1_GeV, 0_ns)};
  const auto forward = response(input);
  std::reverse(input.begin(), input.end());
  const auto reversed = response(input);
  BOOST_REQUIRE_EQUAL(forward.size(), 2);
  BOOST_REQUIRE_EQUAL(reversed.size(), 2);
  BOOST_CHECK_EQUAL(forward[0].cellId, 0);
  BOOST_CHECK_EQUAL(forward[1].cellId, highId);
  BOOST_CHECK_EQUAL(forward[1].energy, reversed[1].energy);
  BOOST_CHECK(response(SimCalorimeterHitContainer{}).empty());
}

BOOST_AUTO_TEST_CASE(RejectMalformedConfiguration) {
  CalorimeterResponse::Config config;
  config.energyScale = 0;
  BOOST_CHECK_THROW(CalorimeterResponse{config}, std::invalid_argument);
  config.energyScale = std::numeric_limits<double>::infinity();
  BOOST_CHECK_THROW(CalorimeterResponse{config}, std::invalid_argument);
  config.energyScale = 1;
  config.energyThreshold = -1_GeV;
  BOOST_CHECK_THROW(CalorimeterResponse{config}, std::invalid_argument);
  config.energyThreshold = 0;
  config.timeMin = 2_ns;
  config.timeMax = 1_ns;
  BOOST_CHECK_THROW(CalorimeterResponse{config}, std::invalid_argument);
  config.timeMin = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_THROW(CalorimeterResponse{config}, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(RejectMalformedDepositsAndOverflow) {
  const CalorimeterResponse response(CalorimeterResponse::Config{});
  auto input = SimCalorimeterHitContainer{deposit(1, -1_GeV, 0_ns)};
  BOOST_CHECK_THROW(response(input), std::invalid_argument);
  input[0].depositedEnergy = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_THROW(response(input), std::invalid_argument);
  input[0] = deposit(1, 1_GeV, std::numeric_limits<double>::infinity());
  BOOST_CHECK_THROW(response(input), std::invalid_argument);
  input[0] = deposit(1, 1_GeV, 0_ns);
  input[0].position.x() = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_THROW(response(input), std::invalid_argument);
  input = {deposit(1, 1_GeV, 0_ns), deposit(1, 1_GeV, 0_ns)};
  input[1].position.y() = 1_mm;
  BOOST_CHECK_THROW(response(input), std::invalid_argument);
  input = {deposit(1, std::numeric_limits<double>::max(), 0_ns),
           deposit(1, std::numeric_limits<double>::max(), 0_ns)};
  BOOST_CHECK_THROW(response(input), std::overflow_error);
  CalorimeterResponse::Config scaled;
  scaled.energyScale = 2;
  input.resize(1);
  BOOST_CHECK_THROW(CalorimeterResponse{scaled}(input), std::overflow_error);
}

BOOST_AUTO_TEST_CASE(AlgorithmUsesTypedEventCollections) {
  using namespace ActsExamples;
  WhiteBoard board;
  AlgorithmContext context(0, 0, board, 0);
  DummySequenceElement source;
  WriteDataHandle<SimCalorimeterHitContainer> write(&source, "SimHits");
  write.initialize("sim_calo_hits");
  write(context, SimCalorimeterHitContainer{deposit(42, 1_GeV, 2_ns)});

  CalorimeterDigitizationAlgorithm::Config config;
  config.inputSimHits = "sim_calo_hits";
  config.outputHits = "calo_hits";
  config.response.energyScale = 2;
  const CalorimeterDigitizationAlgorithm algorithm(config);
  BOOST_CHECK(algorithm.execute(context) == ProcessCode::SUCCESS);

  ReadDataHandle<CalorimeterHitContainer> read(&source, "CaloHits");
  read.initialize(config.outputHits);
  BOOST_REQUIRE_EQUAL(read(context).size(), 1);
  BOOST_CHECK_EQUAL(read(context)[0].energy, 2_GeV);
  BOOST_REQUIRE_EQUAL(read(context)[0].sourceIndices.size(), 1);
  BOOST_CHECK_EQUAL(read(context)[0].sourceIndices[0], 0);
  ReadDataHandle<SimCalorimeterHitContainer> input(&source, "SimHits");
  input.initialize(config.inputSimHits);
  BOOST_CHECK_EQUAL(input(context)[0].depositedEnergy, 1_GeV);

  config.outputHits = config.inputSimHits;
  BOOST_CHECK_THROW(CalorimeterDigitizationAlgorithm{config},
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
