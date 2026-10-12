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
#include "Jets/JetFixture.hpp"

namespace ActsTests {
namespace {
using namespace ActsExamples;
using namespace ActsCalorimeter;
namespace Fixture = ActsExamples::CalorimeterFixture;
void reconstructJets(
    CalorimeterIoEvent& event,
    const CalorimeterJetAlgorithm::Config& config = Fixture::jetConfig()) {
  event.reconstruct();
  event.persistClusters();
  CalorimeterJetAlgorithm{config}.execute(event.context);
}
}  // namespace
BOOST_AUTO_TEST_SUITE(CalorimeterEDM4hepJetsSuite)
BOOST_AUTO_TEST_CASE(JetMomentumMassAndClusterRelations) {
  CalorimeterIoEvent event;
  reconstructJets(event);
  const EDM4hepCalorimeterJetOutputConverter converter{
      Fixture::jetOutputConfig()};
  BOOST_CHECK(converter.execute(event.context) == ProcessCode::SUCCESS);
  const auto& jets =
      event.readPodio<edm4hep::ReconstructedParticleCollection>("CaloJets");
  const auto& clusters =
      event.readPodio<edm4hep::ClusterCollection>("CaloClusters");
  const auto& native = event.read<CalorimeterJetContainer>("calo_jets");
  BOOST_REQUIRE_EQUAL(jets.size(), 2);
  BOOST_REQUIRE_EQUAL(jets[0].clusters_size(), 2);
  BOOST_CHECK(jets[0].getClusters(0) == clusters[0] &&
              jets[0].getClusters(1) == clusters[2]);
  BOOST_CHECK(jets[1].getClusters(0) == clusters[1]);
  BOOST_CHECK_CLOSE(jets[0].getEnergy(), 0.18, 1e-4);
  BOOST_CHECK(jets[0].getMass() > 0);
  BOOST_CHECK_CLOSE(jets[0].getMomentum().x,
                    native[0].fourMomentum[0] / Acts::UnitConstants::GeV, 1e-4);
  const auto metadata = converter.metadata();
  BOOST_CHECK(metadata.numbers.at("acts.calo.jets.radius") ==
              std::vector<double>{0.4});
  BOOST_CHECK(metadata.strings.at("acts.calo.jets.output") ==
              std::vector<std::string>{"CaloJets"});
}
BOOST_AUTO_TEST_CASE(RejectMalformedJetsAndMismatchedConstituents) {
  const std::vector<std::function<void(CalorimeterJetContainer&)>> changes{
      [](auto& jets) { jets[0].clusterIndices.clear(); },
      [](auto& jets) { jets[0].clusterIndices = {999}; },
      [](auto& jets) { jets[0].clusterIndices = {0, 0}; },
      [](auto& jets) { jets[0].clusterIndices = {1}; },
      [](auto& jets) { jets[1].clusterIndices = {0}; },
      [](auto& jets) { jets[0].fourMomentum[0] += 1; },
      [](auto& jets) { jets[0].fourMomentum[3] *= 2; },
      [](auto& jets) { jets[0].fourMomentum[3] = -1; },
      [](auto& jets) {
        jets[0].fourMomentum[2] = std::numeric_limits<double>::quiet_NaN();
      }};
  for (const auto& change : changes) {
    CalorimeterIoEvent event;
    reconstructJets(event);
    auto jets = event.read<CalorimeterJetContainer>("calo_jets");
    change(jets);
    event.write("changed_jets", std::move(jets));
    auto config = Fixture::jetOutputConfig();
    config.inputJets = "changed_jets";
    BOOST_CHECK_THROW(
        EDM4hepCalorimeterJetOutputConverter{config}.execute(event.context),
        std::invalid_argument);
    BOOST_CHECK(!event.board.exists("CaloJets"));
  }
}
BOOST_AUTO_TEST_CASE(RejectJetEnergyOverflowWithRepresentableClusters) {
  CalorimeterIoEvent event;
  event.reconstruct(4e39);
  event.persistClusters();
  CalorimeterJetAlgorithm{Fixture::jetConfig()}.execute(event.context);
  BOOST_CHECK_THROW(
      EDM4hepCalorimeterJetOutputConverter{Fixture::jetOutputConfig()}.execute(
          event.context),
      std::overflow_error);
  BOOST_CHECK(!event.board.exists("CaloJets"));
}

BOOST_AUTO_TEST_CASE(RejectWrongClusterMapping) {
  CalorimeterIoEvent event;
  reconstructJets(event);
  const auto& original =
      event.readPodio<edm4hep::ClusterCollection>("CaloClusters");
  const auto& hits =
      event.readPodio<edm4hep::CalorimeterHitCollection>("CaloHits");
  edm4hep::ClusterCollection changed;
  for (const auto& cluster : original) {
    auto copy = changed.create();
    copy.setEnergy(cluster.getEnergy());
    copy.setPosition(cluster.getPosition());
    for (std::size_t index = 0; index < cluster.hits_size(); ++index) {
      copy.addToHits(hits[0]);
    }
  }
  event.writePodio("wrong_clusters", std::move(changed));
  auto config = Fixture::jetOutputConfig();
  config.inputEdmClusters = "wrong_clusters";
  BOOST_CHECK_THROW(
      EDM4hepCalorimeterJetOutputConverter{config}.execute(event.context),
      std::invalid_argument);
}
BOOST_AUTO_TEST_CASE(DisplacedOriginAndEmptyEvents) {
  CalorimeterIoEvent event;
  auto nativeConfig = Fixture::jetConfig();
  nativeConfig.jets.origin =
      Acts::Vector3{10, 20, 30} * Acts::UnitConstants::mm;
  reconstructJets(event, nativeConfig);
  auto output = Fixture::jetOutputConfig();
  output.reconstruction = nativeConfig.jets;
  const EDM4hepCalorimeterJetOutputConverter converter{output};
  BOOST_CHECK(converter.execute(event.context) == ProcessCode::SUCCESS);
  BOOST_CHECK(converter.metadata().numbers.at("acts.calo.jets.origin") ==
              std::vector<double>({10, 20, 30}));
  CalorimeterIoEvent empty{emptyCalorimeterFrame()};
  reconstructJets(empty);
  EDM4hepCalorimeterJetOutputConverter{Fixture::jetOutputConfig()}.execute(
      empty.context);
  BOOST_CHECK(
      empty.readPodio<edm4hep::ReconstructedParticleCollection>("CaloJets")
          .empty());
}
BOOST_AUTO_TEST_CASE(RejectInvalidJetConfiguration) {
  auto config = Fixture::jetOutputConfig();
  config.outputJets = config.inputJets;
  BOOST_CHECK_THROW(EDM4hepCalorimeterJetOutputConverter{config},
                    std::invalid_argument);
  config = Fixture::jetOutputConfig();
  config.reconstruction.radius = -1;
  BOOST_CHECK_THROW(EDM4hepCalorimeterJetOutputConverter{config},
                    std::invalid_argument);
  CalorimeterIoEvent event;
  reconstructJets(event);
  config = Fixture::jetOutputConfig();
  config.reconstruction.origin.x() = 100 * Acts::UnitConstants::mm;
  BOOST_CHECK_THROW(
      EDM4hepCalorimeterJetOutputConverter{config}.execute(event.context),
      std::invalid_argument);
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
