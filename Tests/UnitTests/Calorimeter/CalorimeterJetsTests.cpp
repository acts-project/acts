// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "ActsCalorimeter/Jets/CalorimeterJetReconstruction.hpp"
#include "ActsExamples/Calorimeter/CalorimeterClusteringAlgorithm.hpp"
#include "ActsExamples/Calorimeter/CalorimeterDigitizationAlgorithm.hpp"
#include "ActsExamples/Calorimeter/CalorimeterJetAlgorithm.hpp"
#include "ActsTests/CommonHelpers/WhiteBoardUtilities.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>

namespace ActsTests {
namespace {
using namespace Acts::UnitLiterals;
using namespace ActsCalorimeter;

CalorimeterCluster cluster(double energy, double phi, double eta = 0) {
  CalorimeterCluster result;
  result.energy = energy;
  result.position =
      Acts::Vector3{1500_mm * std::cos(phi), 1500_mm * std::sin(phi),
                    1500_mm * std::sinh(eta)};
  return result;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(CalorimeterJetsSuite)

BOOST_AUTO_TEST_CASE(AntiKtRecombinationConservesFourMomentumAndLinks) {
  const CalorimeterJetReconstruction reconstruction({});
  const CalorimeterClusterContainer clusters{
      cluster(40_GeV, 0), cluster(0_GeV, 0), cluster(10_GeV, 0.1),
      cluster(15_GeV, 2)};
  const auto jets = reconstruction(clusters);
  BOOST_REQUIRE_EQUAL(jets.size(), 2);
  BOOST_CHECK(jets[0].clusterIndices == std::vector<std::size_t>({0, 2}));
  BOOST_CHECK(jets[1].clusterIndices == std::vector<std::size_t>({3}));
  BOOST_CHECK_CLOSE(jets[0].fourMomentum[3] / 1_GeV, 50, 1e-10);
  const double massSquared = std::pow(jets[0].fourMomentum[3], 2) -
                             jets[0].fourMomentum.head<3>().squaredNorm();
  BOOST_CHECK_CLOSE(massSquared / (1_GeV * 1_GeV),
                    2 * 40 * 10 * (1 - std::cos(0.1)), 1e-8);
  Acts::Vector4 sum = Acts::Vector4::Zero();
  std::set<std::size_t> used;
  for (const auto& jet : jets) {
    sum += jet.fourMomentum;
    for (const auto index : jet.clusterIndices) {
      BOOST_CHECK(used.insert(index).second);
    }
  }
  const Acts::Vector4 expected{
      (40 + 10 * std::cos(0.1) + 15 * std::cos(2)) * 1_GeV,
      (10 * std::sin(0.1) + 15 * std::sin(2)) * 1_GeV, 0, 65_GeV};
  BOOST_CHECK_SMALL((sum - expected).norm() / 1_GeV, 1e-12);
  BOOST_CHECK_EQUAL(used.size(), 3);
  BOOST_CHECK(!used.contains(1));
  // Returned values remain valid after the local FastJet sequence is gone.
  BOOST_CHECK_EQUAL(jets[1].fourMomentum[3], 15_GeV);
  BOOST_CHECK_EQUAL(clusters[0].energy, 40_GeV);
}

BOOST_AUTO_TEST_CASE(PhiWrapAndConfigurableRadius) {
  const CalorimeterClusterContainer clusters{
      cluster(10_GeV, std::numbers::pi - 0.05),
      cluster(5_GeV, -std::numbers::pi + 0.05)};
  const auto joined = CalorimeterJetReconstruction{{0.4, 0}}(clusters);
  BOOST_REQUIRE_EQUAL(joined.size(), 1);
  BOOST_CHECK(joined[0].clusterIndices == std::vector<std::size_t>({0, 1}));
  BOOST_CHECK_EQUAL(joined[0].fourMomentum[3], 15_GeV);
  const auto split = CalorimeterJetReconstruction{{0.05, 0}}(clusters);
  BOOST_REQUIRE_EQUAL(split.size(), 2);
  BOOST_CHECK(split[0].clusterIndices == std::vector<std::size_t>({0}));
}

BOOST_AUTO_TEST_CASE(NonzeroRapidityAndPtOrdering) {
  const CalorimeterClusterContainer clusters{cluster(20_GeV, 0, 2),
                                             cluster(10_GeV, 2, 0)};
  const auto jets = CalorimeterJetReconstruction{{}}(clusters);
  BOOST_REQUIRE_EQUAL(jets.size(), 2);
  // Energy ordering differs from pT ordering for forward clusters.
  BOOST_CHECK(jets[0].clusterIndices == std::vector<std::size_t>({1}));
  BOOST_CHECK(jets[1].clusterIndices == std::vector<std::size_t>({0}));
  BOOST_CHECK_CLOSE(jets[1].fourMomentum[0] / 1_GeV, 20 / std::cosh(2), 1e-10);
  BOOST_CHECK_CLOSE(jets[1].fourMomentum[2] / 1_GeV, 20 * std::tanh(2), 1e-10);
}

BOOST_AUTO_TEST_CASE(DisplacedOriginAndInclusivePtCut) {
  CalorimeterJetReconstruction::Config config;
  config.origin = Acts::Vector3{10_mm, 20_mm, 30_mm};
  config.jetPtMin = 3_GeV;
  auto input = cluster(5_GeV, 0);
  input.position = config.origin + Acts::Vector3{3_mm, 0, 4_mm};
  const CalorimeterClusterContainer clusters{input};
  const auto jets = CalorimeterJetReconstruction{config}(clusters);
  BOOST_REQUIRE_EQUAL(jets.size(), 1);
  BOOST_CHECK_SMALL(
      (jets[0].fourMomentum - Acts::Vector4{3_GeV, 0, 4_GeV, 5_GeV}).norm() /
          1_GeV,
      1e-12);
  config.jetPtMin =
      std::nextafter(3_GeV, std::numeric_limits<double>::infinity());
  BOOST_CHECK(CalorimeterJetReconstruction{config}(clusters).empty());
  config.jetPtMin = std::numeric_limits<double>::max();
  BOOST_CHECK(CalorimeterJetReconstruction{config}(clusters).empty());
}

BOOST_AUTO_TEST_CASE(EmptyEventsZeroEnergyAndDeterministicTies) {
  const CalorimeterJetReconstruction reconstruction({});
  BOOST_CHECK(reconstruction(CalorimeterClusterContainer{}).empty());
  BOOST_CHECK(reconstruction(CalorimeterClusterContainer{{}}).empty());
  auto negativeX = cluster(10_GeV, 0);
  negativeX.position.x() = -1500_mm;
  const CalorimeterClusterContainer clusters{cluster(10_GeV, 0), negativeX};
  const auto jets = reconstruction(clusters);
  BOOST_REQUIRE_EQUAL(jets.size(), 2);
  BOOST_CHECK(jets[0].clusterIndices == std::vector<std::size_t>({0}));
  BOOST_CHECK(jets[1].clusterIndices == std::vector<std::size_t>({1}));
  BOOST_CHECK(reconstruction(CalorimeterClusterContainer{}).empty());
  BOOST_CHECK_EQUAL(reconstruction(clusters).size(), 2);
}

BOOST_AUTO_TEST_CASE(RejectInvalidConfiguration) {
  CalorimeterJetReconstruction::Config config;
  for (double value : {0.0, -1.0, 1e6, std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity()}) {
    config.radius = value;
    BOOST_CHECK_THROW(CalorimeterJetReconstruction{config},
                      std::invalid_argument);
  }
  config.radius = 0.4;
  for (double value : {-1.0, std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity()}) {
    config.jetPtMin = value;
    BOOST_CHECK_THROW(CalorimeterJetReconstruction{config},
                      std::invalid_argument);
  }
  config.jetPtMin = 0;
  config.origin.x() = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_THROW(CalorimeterJetReconstruction{config},
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(RejectMalformedKinematicsAndNumericOverflow) {
  const CalorimeterJetReconstruction reconstruction({});
  CalorimeterClusterContainer clusters{cluster(1_GeV, 0)};
  for (double value : {-1_GeV, std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity()}) {
    clusters[0].energy = value;
    BOOST_CHECK_THROW(reconstruction(clusters), std::invalid_argument);
  }
  clusters[0] = cluster(1_GeV, 0);
  clusters[0].position.x() = std::numeric_limits<double>::infinity();
  BOOST_CHECK_THROW(reconstruction(clusters), std::invalid_argument);
  clusters[0].energy = 0;
  BOOST_CHECK_THROW(reconstruction(clusters), std::invalid_argument);
  clusters[0] = cluster(1_GeV, 0);
  clusters[0].position = Acts::Vector3::Zero();
  BOOST_CHECK_THROW(reconstruction(clusters), std::invalid_argument);
  clusters[0] = cluster(std::numeric_limits<double>::max(), 0);
  BOOST_CHECK_THROW(reconstruction(clusters), std::overflow_error);
  const double largeEnergy = std::sqrt(std::numeric_limits<double>::max()) / 3;
  clusters = {cluster(largeEnergy, 0), cluster(largeEnergy, 0.1)};
  BOOST_CHECK_THROW(reconstruction(clusters), std::overflow_error);
}

BOOST_AUTO_TEST_CASE(
    AlgorithmPreservesProvenanceAcrossResponseClustersAndJets) {
  using namespace ActsExamples;
  WhiteBoard board;
  AlgorithmContext context(0, 0, board, 0);
  DummySequenceElement source;
  WriteDataHandle<SimCalorimeterHitContainer> write(&source, "Deposits");
  write.initialize("sim_calo_hits");
  write(context, SimCalorimeterHitContainer{
                     {42, Acts::Vector3{1500_mm, 0, 0}, 5_GeV, 1_ns},
                     {42, Acts::Vector3{1500_mm, 0, 0}, 2_GeV, 2_ns},
                     {43, Acts::Vector3{1500_mm, 10_mm, 0}, 3_GeV, 1_ns}});
  CalorimeterDigitizationAlgorithm::Config response;
  response.inputSimHits = "sim_calo_hits";
  response.outputHits = "calo_hits";
  response.response.energyScale = 2;
  BOOST_CHECK(CalorimeterDigitizationAlgorithm{response}.execute(context) ==
              ProcessCode::SUCCESS);
  CalorimeterClusteringAlgorithm::Config clustering;
  clustering.inputHits = "calo_hits";
  clustering.outputClusters = "calo_clusters";
  // No neighbours: two clusters later merge into one jet.
  BOOST_CHECK(CalorimeterClusteringAlgorithm{clustering}.execute(context) ==
              ProcessCode::SUCCESS);
  CalorimeterJetAlgorithm::Config config;
  config.inputClusters = "calo_clusters";
  config.outputJets = "calo_jets";
  const CalorimeterJetAlgorithm algorithm(config);
  BOOST_CHECK(algorithm.execute(context) == ProcessCode::SUCCESS);
  ReadDataHandle<CalorimeterJetContainer> read(&source, "Jets");
  read.initialize(config.outputJets);
  ReadDataHandle<CalorimeterClusterContainer> clusters(&source, "Clusters");
  clusters.initialize(config.inputClusters);
  ReadDataHandle<CalorimeterHitContainer> hits(&source, "Hits");
  hits.initialize(clustering.inputHits);
  ReadDataHandle<SimCalorimeterHitContainer> deposits(&source, "Deposits");
  deposits.initialize(response.inputSimHits);
  const auto& jets = read(context);
  BOOST_REQUIRE_EQUAL(jets.size(), 1);
  BOOST_CHECK_EQUAL(jets[0].fourMomentum[3], 20_GeV);
  BOOST_CHECK(jets[0].clusterIndices == std::vector<std::size_t>({0, 1}));
  std::set<std::size_t> reached;
  double depositedEnergy = 0;
  for (const auto clusterIndex : jets[0].clusterIndices) {
    for (const auto hitIndex : clusters(context).at(clusterIndex).hitIndices) {
      for (const auto depositIndex : hits(context).at(hitIndex).sourceIndices) {
        BOOST_CHECK(reached.insert(depositIndex).second);
        depositedEnergy += deposits(context).at(depositIndex).depositedEnergy;
      }
    }
  }
  BOOST_CHECK_EQUAL(reached.size(), 3);
  BOOST_CHECK_EQUAL(depositedEnergy, 10_GeV);
  BOOST_CHECK_EQUAL(hits(context)[0].energy, 14_GeV);
  BOOST_CHECK_EQUAL(clusters(context)[0].energy, 14_GeV);
  config.outputJets = config.inputClusters;
  BOOST_CHECK_THROW(CalorimeterJetAlgorithm{config}, std::invalid_argument);
  config.outputJets.clear();
  BOOST_CHECK_THROW(CalorimeterJetAlgorithm{config}, std::invalid_argument);
  config.outputJets = "jets";
  config.inputClusters.clear();
  BOOST_CHECK_THROW(CalorimeterJetAlgorithm{config}, std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
