// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "ActsCalorimeter/Clustering/CalorimeterClusterer.hpp"
#include "ActsExamples/Calorimeter/CalorimeterClusteringAlgorithm.hpp"
#include "ActsTests/CommonHelpers/WhiteBoardUtilities.hpp"

#include <algorithm>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

namespace ActsTests {
namespace {
using namespace Acts::UnitLiterals;
using namespace ActsCalorimeter;

CalorimeterHit hit(std::uint64_t id, double energy, double y = 0_mm,
                   double time = 0_ns) {
  return {id, Acts::Vector3{1500_mm, y, 0}, energy, time, {}};
}

std::vector<std::uint64_t> constituentIds(const CalorimeterCluster& cluster,
                                          const CalorimeterHitContainer& hits) {
  std::vector<std::uint64_t> ids;
  for (const auto index : cluster.hitIndices) {
    ids.push_back(hits.at(index).cellId);
  }
  return ids;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(CalorimeterClusteringSuite)

BOOST_AUTO_TEST_CASE(SeededComponentsAndEnergyAccounting) {
  const CellNeighbourContainer edges{{1, 2}, {2, 3}, {4, 5}};
  const CalorimeterClusterer clusterer({1_GeV, 0.1_GeV}, edges);
  const CalorimeterHitContainer hits{hit(1, 2_GeV, 0_mm, 1_ns),
                                     hit(2, 0.4_GeV, 10_mm, 4_ns),
                                     hit(3, 0.1_GeV, 20_mm, 6_ns),
                                     hit(4, 0.8_GeV),
                                     hit(5, 0.8_GeV),
                                     hit(6, 1_GeV)};
  const auto clusters = clusterer(hits);
  BOOST_REQUIRE_EQUAL(clusters.size(), 2);
  BOOST_CHECK_CLOSE(clusters[0].energy / 1_GeV, 2.5, 1e-10);
  BOOST_CHECK_CLOSE(clusters[0].position.y() / 1_mm, 2.4, 1e-10);
  BOOST_CHECK_CLOSE(clusters[0].time / 1_ns, 1.68, 1e-10);
  BOOST_CHECK_EQUAL(clusters[0].seedCellId, 1);
  BOOST_CHECK_EQUAL(clusters[1].seedCellId, 6);
  BOOST_CHECK_EQUAL(clusters[1].energy, 1_GeV);

  double acceptedEnergy = 0;
  std::set<std::size_t> used;
  for (const auto& cluster : clusters) {
    double constituentEnergy = 0;
    for (const auto index : cluster.hitIndices) {
      BOOST_CHECK(used.insert(index).second);
      constituentEnergy += hits.at(index).energy;
    }
    BOOST_CHECK_EQUAL(cluster.energy, constituentEnergy);
    acceptedEnergy += cluster.energy;
  }
  BOOST_CHECK_CLOSE(acceptedEnergy / 1_GeV, 3.5, 1e-10);
  BOOST_CHECK_EQUAL(used.size(), 4);
  // The seed cut is on an individual cell, not the component's energy sum.
  BOOST_CHECK(!used.contains(3));
  BOOST_CHECK(!used.contains(4));
}

BOOST_AUTO_TEST_CASE(BelowThresholdOrMissingCellsDoNotBridge) {
  const CellNeighbourContainer edges{{1, 2}, {2, 3}};
  const CalorimeterClusterer clusterer({1_GeV, 0.1_GeV}, edges);
  CalorimeterHitContainer hits{hit(1, 2_GeV), hit(2, 0.05_GeV), hit(3, 2_GeV)};
  BOOST_CHECK_EQUAL(clusterer(hits).size(), 2);
  hits[1].energy = 0;
  BOOST_CHECK_EQUAL(clusterer(hits).size(), 2);
  hits.erase(hits.begin() + 1);
  BOOST_CHECK_EQUAL(clusterer(hits).size(), 2);
}

BOOST_AUTO_TEST_CASE(MultipleSeedsAndBoundaryCellsCountOnce) {
  const CellNeighbourContainer edges{{1, 2}, {2, 3}, {3, 4}, {4, 1}};
  const CalorimeterClusterer clusterer({0.5_GeV, 0.1_GeV}, edges);
  const CalorimeterHitContainer hits{hit(1, 0.5_GeV), hit(2, 0.1_GeV),
                                     hit(3, 0.5_GeV), hit(4, 0_GeV)};
  const auto clusters = clusterer(hits);
  BOOST_REQUIRE_EQUAL(clusters.size(), 1);
  BOOST_CHECK_CLOSE(clusters[0].energy / 1_GeV, 1.1, 1e-10);
  BOOST_CHECK_EQUAL(clusters[0].seedCellId, 1);
  BOOST_CHECK_EQUAL(clusters[0].hitIndices.size(), 3);
}

BOOST_AUTO_TEST_CASE(TopologyAndInputPermutationPreserveResults) {
  const auto highestId = std::numeric_limits<std::uint64_t>::max();
  CellNeighbourContainer edges{
      {highestId, 0}, {0, highestId}, {highestId, 0}, {10, 11}};
  CalorimeterHitContainer hits{hit(highestId, 1_GeV, 2_mm, 3_ns),
                               hit(0, 2_GeV, 4_mm, 5_ns), hit(7, 1_GeV)};
  const CalorimeterClusterer forward({}, edges);
  const auto original = forward(hits);
  std::reverse(edges.begin(), edges.end());
  auto permutedHits = hits;
  std::reverse(permutedHits.begin(), permutedHits.end());
  const CalorimeterClusterer reverse({}, edges);
  const auto permuted = reverse(permutedHits);
  BOOST_REQUIRE_EQUAL(original.size(), 2);
  BOOST_REQUIRE_EQUAL(permuted.size(), 2);
  for (std::size_t i = 0; i < original.size(); ++i) {
    BOOST_CHECK_EQUAL(original[i].energy, permuted[i].energy);
    BOOST_CHECK_EQUAL(original[i].time, permuted[i].time);
    BOOST_CHECK(original[i].position == permuted[i].position);
    BOOST_CHECK_EQUAL(original[i].seedCellId, permuted[i].seedCellId);
    BOOST_CHECK(constituentIds(original[i], hits) ==
                constituentIds(permuted[i], permutedHits));
  }
  BOOST_CHECK(constituentIds(original[0], hits) ==
              std::vector<std::uint64_t>({0, highestId}));
  // Numeric closeness of IDs has no adjacency meaning.
  BOOST_CHECK_EQUAL(CalorimeterClusterer{{}}(hits).size(), 3);
}

BOOST_AUTO_TEST_CASE(EmptyEventsAndImmutableTopology) {
  CellNeighbourContainer edges{{1, 2}};
  const CalorimeterClusterer clusterer({}, edges);
  edges.clear();
  const CalorimeterHitContainer hits{hit(1, 1_GeV), hit(2, 1_GeV)};
  BOOST_CHECK_EQUAL(clusterer(hits).size(), 1);
  BOOST_CHECK(clusterer(CalorimeterHitContainer{}).empty());
  BOOST_CHECK_EQUAL(clusterer(hits).size(), 1);
  BOOST_CHECK(
      CalorimeterClusterer{{}}(CalorimeterHitContainer{hit(1, 0_GeV)}).empty());
}

BOOST_AUTO_TEST_CASE(RejectInvalidConfigurationAndTopology) {
  BOOST_CHECK_THROW(CalorimeterClusterer({0_GeV, 1_GeV}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(CalorimeterClusterer({1_GeV, -1_GeV}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(
      CalorimeterClusterer({std::numeric_limits<double>::quiet_NaN(), 0}),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      CalorimeterClusterer({std::numeric_limits<double>::infinity(), 0}),
      std::invalid_argument);
  const CellNeighbourContainer selfEdge{{42, 42}};
  BOOST_CHECK_THROW(CalorimeterClusterer({}, selfEdge), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(RejectMalformedCellsAndOverflow) {
  const CellNeighbourContainer edges{{1, 2}};
  const CalorimeterClusterer clusterer({}, edges);
  auto hits = CalorimeterHitContainer{hit(1, 1_GeV), hit(1, 0_GeV)};
  BOOST_CHECK_THROW(clusterer(hits), std::invalid_argument);
  hits = {hit(1, -1_GeV)};
  BOOST_CHECK_THROW(clusterer(hits), std::invalid_argument);
  hits[0].energy = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_THROW(clusterer(hits), std::invalid_argument);
  hits[0] = hit(1, 1_GeV);
  hits[0].time = std::numeric_limits<double>::infinity();
  BOOST_CHECK_THROW(clusterer(hits), std::invalid_argument);
  hits[0] = hit(1, 1_GeV);
  hits[0].position.x() = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_THROW(clusterer(hits), std::invalid_argument);
  hits = {hit(1, std::numeric_limits<double>::max()),
          hit(2, std::numeric_limits<double>::max())};
  BOOST_CHECK_THROW(clusterer(hits), std::overflow_error);
}

BOOST_AUTO_TEST_CASE(AlgorithmPreservesHitAndDepositProvenance) {
  using namespace ActsExamples;
  WhiteBoard board;
  AlgorithmContext context(0, 0, board, 0);
  DummySequenceElement source;
  WriteDataHandle<CalorimeterHitContainer> write(&source, "Hits");
  write.initialize("calo_hits");
  auto cells = CalorimeterHitContainer{hit(2, 1_GeV), hit(1, 2_GeV)};
  cells[0].sourceIndices = {4};
  cells[1].sourceIndices = {0, 3};
  write(context, std::move(cells));

  CalorimeterClusteringAlgorithm::Config config;
  config.inputHits = "calo_hits";
  config.outputClusters = "calo_clusters";
  config.neighbours = {{1, 2}};
  const CalorimeterClusteringAlgorithm algorithm(config);
  BOOST_CHECK(algorithm.execute(context) == ProcessCode::SUCCESS);
  ReadDataHandle<CalorimeterClusterContainer> read(&source, "Clusters");
  read.initialize(config.outputClusters);
  ReadDataHandle<CalorimeterHitContainer> input(&source, "Hits");
  input.initialize(config.inputHits);
  const auto& clusters = read(context);
  BOOST_REQUIRE_EQUAL(clusters.size(), 1);
  BOOST_CHECK_EQUAL(clusters[0].energy, 3_GeV);
  BOOST_REQUIRE_EQUAL(clusters[0].hitIndices.size(), 2);
  BOOST_CHECK_EQUAL(clusters[0].hitIndices[0], 1);
  BOOST_CHECK(input(context)[clusters[0].hitIndices[0]].sourceIndices ==
              std::vector<std::size_t>({0, 3}));
  BOOST_CHECK_EQUAL(input(context)[0].energy, 1_GeV);
  config.outputClusters = config.inputHits;
  BOOST_CHECK_THROW(CalorimeterClusteringAlgorithm{config},
                    std::invalid_argument);
  config.inputHits.clear();
  BOOST_CHECK_THROW(CalorimeterClusteringAlgorithm{config},
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
