// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

#include "EDM4hepReconstructionTestHelpers.hpp"

namespace ActsTests {
namespace {
using namespace ActsExamples;
using namespace ActsCalorimeter;
namespace Fixture = ActsExamples::CalorimeterFixture;
}  // namespace
BOOST_AUTO_TEST_SUITE(CalorimeterEDM4hepReconstructionSuite)
BOOST_AUTO_TEST_CASE(ClusterRelationsTimesAndSeedIds) {
  CalorimeterIoEvent event;
  event.reconstruct();
  event.persistClusters();
  const auto& clusters =
      event.readPodio<edm4hep::ClusterCollection>("CaloClusters");
  const auto& times =
      event.readPodio<podio::UserDataCollection<double>>("CaloClusterTimes");
  const auto& seeds = event.readPodio<podio::UserDataCollection<std::uint64_t>>(
      "CaloClusterSeedCellIds");
  const auto& hits =
      event.readPodio<edm4hep::CalorimeterHitCollection>("CaloHits");
  BOOST_REQUIRE_EQUAL(clusters.size(), 3);
  BOOST_REQUIRE_EQUAL(times.size(), 3);
  BOOST_REQUIRE_EQUAL(seeds.size(), 3);
  BOOST_CHECK_CLOSE(clusters[2].getEnergy(), 0.12, 1e-4);
  BOOST_REQUIRE_EQUAL(clusters[2].hits_size(), 2);
  BOOST_CHECK(clusters[2].getHits(0) == hits[2]);
  BOOST_CHECK(clusters[2].getHits(1) == hits[3]);
  BOOST_CHECK_EQUAL(seeds[2], Fixture::highCellId);
  BOOST_CHECK_CLOSE(times[2], 10. / 3, 1e-4);
  BOOST_CHECK_CLOSE(clusters[2].getPosition().y, 40. / 3, 1e-4);
  BOOST_CHECK(event.board.exists("calo_hits") &&
              event.board.exists("calo_clusters"));
}
BOOST_AUTO_TEST_CASE(RejectMalformedClustersWithoutWritingPartialOutput) {
  const std::vector<std::function<void(CalorimeterClusterContainer&)>> changes{
      [](auto& clusters) { clusters[2].hitIndices.clear(); },
      [](auto& clusters) { clusters[2].hitIndices = {999}; },
      [](auto& clusters) { clusters[2].hitIndices = {2, 2}; },
      [](auto& clusters) { clusters[2].hitIndices = {0}; },
      [](auto& clusters) { clusters[2].energy *= 2; },
      [](auto& clusters) { clusters[2].seedCellId = 7; },
      [](auto& clusters) { clusters[2].time += 1; },
      [](auto& clusters) { clusters[2].position.x() += 1; },
      [](auto& clusters) {
        clusters[2].energy = std::numeric_limits<double>::infinity();
      },
      [](auto& clusters) {
        clusters[2].position.y() = std::numeric_limits<double>::quiet_NaN();
      }};
  for (const auto& change : changes) {
    CalorimeterIoEvent event;
    event.reconstruct();
    auto native = event.read<CalorimeterClusterContainer>("calo_clusters");
    change(native);
    event.write("changed_clusters", std::move(native));
    auto config = Fixture::clusterOutputConfig();
    config.inputClusters = "changed_clusters";
    BOOST_CHECK_THROW(
        EDM4hepCalorimeterClusterOutputConverter{config}.execute(event.context),
        std::invalid_argument);
    BOOST_CHECK(!event.board.exists("CaloClusters") &&
                !event.board.exists("CaloClusterTimes") &&
                !event.board.exists("CaloClusterSeedCellIds"));
  }
}
BOOST_AUTO_TEST_CASE(RejectClusterEnergyOverflowWithRepresentableHits) {
  CalorimeterIoEvent event;
  event.reconstruct(6e39);
  BOOST_CHECK_THROW(event.persistClusters(), std::overflow_error);
  BOOST_CHECK(!event.board.exists("CaloClusters"));
}

BOOST_AUTO_TEST_CASE(RejectWrongPersistentHitMapping) {
  for (bool empty : {false, true}) {
    CalorimeterIoEvent event;
    event.reconstruct();
    edm4hep::CalorimeterHitCollection changed;
    if (!empty) {
      const auto& original =
          event.readPodio<edm4hep::CalorimeterHitCollection>("CaloHits");
      for (std::size_t index = 0; index < original.size(); ++index) {
        const auto hit = original[(index + 1) % original.size()];
        auto copy = changed.create();
        copy.setCellID(hit.getCellID());
        copy.setEnergy(hit.getEnergy());
        copy.setTime(hit.getTime());
        copy.setPosition(hit.getPosition());
      }
    }
    event.writePodio("wrong_hits", std::move(changed));
    auto config = Fixture::clusterOutputConfig();
    config.inputEdmHits = "wrong_hits";
    BOOST_CHECK_THROW(
        EDM4hepCalorimeterClusterOutputConverter{config}.execute(event.context),
        std::invalid_argument);
  }
}
BOOST_AUTO_TEST_CASE(EmptyClusterCollectionsAndMetadata) {
  CalorimeterIoEvent event{emptyCalorimeterFrame()};
  event.reconstruct();
  event.persistClusters();
  BOOST_CHECK(
      event.readPodio<edm4hep::ClusterCollection>("CaloClusters").empty());
  BOOST_CHECK(
      event.readPodio<podio::UserDataCollection<double>>("CaloClusterTimes")
          .empty());
  EDM4hepCalorimeterMetadata{Fixture::metadataConfig()}.execute(event.context);
  const auto& frame = event.read<podio::Frame>("calo_events");
  BOOST_CHECK(frame.getParameter<double>("acts.calo.schemaVersion") == 1);
  BOOST_CHECK(!event.board.exists("events"));
}
BOOST_AUTO_TEST_CASE(MetadataPreservesFrameAndRecordsActualConfiguration) {
  auto frame = Fixture::makeFrame(true);
  frame.putParameter("generator.tag", std::string("existing value"));
  CalorimeterIoEvent event{std::move(frame)};
  event.reconstruct();
  event.persistClusters();
  auto config = Fixture::metadataConfig();
  config.clustering.neighbours.push_back(
      {Fixture::highCellId + 1, Fixture::highCellId});
  config.clustering.neighbours.push_back(
      {Fixture::highCellId, Fixture::highCellId + 1});
  BOOST_CHECK(EDM4hepCalorimeterMetadata{config}.execute(event.context) ==
              ProcessCode::SUCCESS);
  const auto& output = event.read<podio::Frame>("calo_events");
  BOOST_CHECK(output.getParameter<std::string>("generator.tag") ==
              "existing value");
  BOOST_CHECK_EQUAL(
      output.get<edm4hep::SimCalorimeterHitCollection>("SimCaloHits").size(),
      5);
  BOOST_CHECK(output.getParameter<double>("acts.calo.response.energyScale") ==
              2);
  BOOST_CHECK(output.getParameter<std::vector<double>>(
                  "acts.calo.response.timeWindow") ==
              std::vector<double>({0, 10}));
  BOOST_CHECK(
      output.getParameter<std::vector<std::string>>(
          "acts.calo.clustering.neighbours") ==
      std::vector<std::string>{std::to_string(Fixture::highCellId) + ":" +
                               std::to_string(Fixture::highCellId + 1)});
  BOOST_CHECK(
      output.getParameter<std::string>("acts.calo.output.clusterSeedCellIds") ==
      "CaloClusterSeedCellIds");
}
BOOST_AUTO_TEST_CASE(RejectMetadataCollisionsAcrossParameterTypes) {
  for (int type = 0; type < 4; ++type) {
    auto frame = Fixture::makeFrame();
    const std::string key = "acts.calo.old";
    if (type == 0) {
      frame.putParameter(key, 1);
    }
    if (type == 1) {
      frame.putParameter(key, 1.f);
    }
    if (type == 2) {
      frame.putParameter(key, 1.);
    }
    if (type == 3) {
      frame.putParameter(key, std::string("old"));
    }
    CalorimeterIoEvent event{std::move(frame)};
    BOOST_CHECK_THROW(
        EDM4hepCalorimeterMetadata{Fixture::metadataConfig()}.execute(
            event.context),
        std::invalid_argument);
    BOOST_CHECK(!event.board.exists("calo_events"));
  }
}
BOOST_AUTO_TEST_CASE(RejectInvalidOutputAndMetadataConfiguration) {
  auto cluster = Fixture::clusterOutputConfig();
  cluster.outputTimes = cluster.outputClusters;
  BOOST_CHECK_THROW(EDM4hepCalorimeterClusterOutputConverter{cluster},
                    std::invalid_argument);
  cluster = Fixture::clusterOutputConfig();
  cluster.inputEdmHits.clear();
  BOOST_CHECK_THROW(EDM4hepCalorimeterClusterOutputConverter{cluster},
                    std::invalid_argument);
  const std::vector<std::function<void(EDM4hepCalorimeterMetadata::Config&)>>
      changes{[](auto& config) { config.outputFrame = config.inputFrame; },
              [](auto& config) { config.geometryIdentifier.clear(); },
              [](auto& config) { config.clustering.inputHits = "wrong_hits"; },
              [](auto& config) {
                config.additional.strings["acts.calo.jets.inputClusters"] = {
                    "wrong_clusters"};
              },
              [](auto& config) {
                config.clustering.clustering.seedEnergyThreshold = -1;
              },
              [](auto& config) {
                config.additional.numbers["acts.calo.schemaVersion"] = {2};
              },
              [](auto& config) {
                config.additional.strings["other.key"] = {"value"};
              },
              [](auto& config) {
                config.additional.numbers["acts.calo.extra"] = {
                    std::numeric_limits<double>::infinity()};
              },
              [](auto& config) {
                config.additional.strings["acts.calo.extra"] = {"value"};
                config.additional.numbers["acts.calo.extra"] = {1};
              }};
  for (const auto& change : changes) {
    auto config = Fixture::metadataConfig();
    change(config);
    BOOST_CHECK_THROW(EDM4hepCalorimeterMetadata{config},
                      std::invalid_argument);
  }
  auto unbounded = Fixture::metadataConfig();
  unbounded.response.response.timeMin =
      -std::numeric_limits<double>::infinity();
  unbounded.response.response.timeMax = std::numeric_limits<double>::infinity();
  CalorimeterIoEvent event;
  EDM4hepCalorimeterMetadata{unbounded}.execute(event.context);
  const auto bounds =
      event.read<podio::Frame>("calo_events")
          .getParameter<std::vector<double>>("acts.calo.response.timeWindow");
  BOOST_REQUIRE(bounds.has_value());
  BOOST_CHECK(std::isinf((*bounds)[0]) && (*bounds)[0] < 0 &&
              std::isinf((*bounds)[1]) && (*bounds)[1] > 0);
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
