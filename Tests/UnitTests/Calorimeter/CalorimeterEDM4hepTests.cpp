// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Calorimeter/CalorimeterDigitizationAlgorithm.hpp"
#include "ActsExamples/Calorimeter/EDM4hepCalorimeterInputConverter.hpp"
#include "ActsExamples/Calorimeter/EDM4hepCalorimeterOutputConverter.hpp"
#include "ActsTests/CommonHelpers/WhiteBoardUtilities.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

#include "ExampleFixture.hpp"

namespace ActsTests {
namespace {
using namespace Acts::UnitLiterals;
using namespace ActsExamples;
using namespace ActsCalorimeter;

EDM4hepCalorimeterInputConverter::Config inputConfig() {
  EDM4hepCalorimeterInputConverter::Config config;
  config.inputSimHits = "SimCaloHits";
  config.outputDeposits = "calo_deposits";
  config.outputSources = "calo_sources";
  config.cellCentre = CalorimeterFixture::cellCentre;
  return config;
}

EDM4hepCalorimeterOutputConverter::Config outputConfig() {
  return {"calo_hits", "calo_deposits", "calo_sources", "CaloHits",
          "CaloLinks"};
}

struct Event {
  WhiteBoard board;
  AlgorithmContext context{0, 0, board, 0};
  DummySequenceElement source;

  explicit Event(podio::Frame frame = CalorimeterFixture::makeFrame()) {
    WriteDataHandle<podio::Frame> write(&source, "Frame");
    write.initialize("events");
    write(context, std::move(frame));
  }

  template <typename T>
  const T& read(const std::string& name) {
    ReadDataHandle<T> handle(&source, "Read");
    handle.initialize(name);
    return handle(context);
  }

  template <typename T>
  void write(const std::string& name, T data) {
    WriteDataHandle<T> handle(&source, "Write");
    handle.initialize(name);
    handle(context, std::move(data));
  }

  void convert() {
    BOOST_CHECK(EDM4hepCalorimeterInputConverter{inputConfig()}.execute(
                    context) == ProcessCode::SUCCESS);
    CalorimeterDigitizationAlgorithm::Config response;
    response.inputSimHits = "calo_deposits";
    response.outputHits = "calo_hits";
    response.response.energyScale = 2;
    response.response.timeMin = 0_ns;
    response.response.timeMax = 10_ns;
    BOOST_CHECK(CalorimeterDigitizationAlgorithm{response}.execute(context) ==
                ProcessCode::SUCCESS);
  }

  template <typename T>
  const T& readPodio(const std::string& name) {
    PodioCollectionReadHandle<T> handle(&source, "ReadPodio");
    handle.initialize(name);
    return handle(context);
  }
};

podio::Frame oneHit(float energy, std::vector<std::pair<float, float>> values,
                    bool repeat = false) {
  edm4hep::SimCalorimeterHitCollection hits;
  edm4hep::CaloHitContributionCollection contributions;
  auto hit = hits.create();
  hit.setCellID(42);
  hit.setEnergy(energy);
  for (const auto [depositEnergy, time] : values) {
    auto contribution = contributions.create();
    contribution.setEnergy(depositEnergy);
    contribution.setTime(time);
    hit.addToContributions(contribution);
    if (repeat) {
      hit.addToContributions(contribution);
    }
  }
  podio::Frame frame;
  frame.put(std::move(contributions), "CaloContributions");
  frame.put(std::move(hits), "SimCaloHits");
  return frame;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(CalorimeterEDM4hepSuite)

BOOST_AUTO_TEST_CASE(ExpandContributionsWithUnitsAndGeometryCentres) {
  Event event;
  auto config = inputConfig();
  std::map<std::uint64_t, int> calls;
  config.cellCentre = [&calls](std::uint64_t id) {
    ++calls[id];
    return CalorimeterFixture::cellCentre(id);
  };
  BOOST_CHECK(EDM4hepCalorimeterInputConverter{config}.execute(event.context) ==
              ProcessCode::SUCCESS);
  const auto& deposits =
      event.read<SimCalorimeterHitContainer>(config.outputDeposits);
  const auto& sources =
      event.read<EDM4hepCalorimeterSourceContainer>(config.outputSources);
  BOOST_REQUIRE_EQUAL(deposits.size(), 5);
  BOOST_REQUIRE_EQUAL(sources.size(), 5);
  BOOST_CHECK_EQUAL(deposits[0].cellId, CalorimeterFixture::highCellId);
  BOOST_CHECK(deposits[0].position ==
              CalorimeterFixture::cellCentre(deposits[0].cellId));
  BOOST_CHECK_CLOSE(deposits[0].depositedEnergy / 1_GeV, 0.01, 1e-4);
  BOOST_CHECK_EQUAL(deposits[2].time, 100_ns);
  BOOST_CHECK_EQUAL(deposits[3].time, 7_ns);
  BOOST_CHECK_EQUAL(calls.at(CalorimeterFixture::highCellId), 1);
  BOOST_CHECK_EQUAL(calls.at(42), 1);
  BOOST_CHECK(sources[0].hit == sources[2].hit);
  BOOST_CHECK(sources[0].hit != sources[3].hit);
  BOOST_CHECK_EQUAL(sources[0].contribution.getParticle().getPDG(), 211);
}

BOOST_AUTO_TEST_CASE(PersistentTruthLinksUseAcceptedEnergyFractions) {
  Event event;
  event.convert();
  const auto config = outputConfig();
  const EDM4hepCalorimeterOutputConverter converter(config);
  BOOST_CHECK(converter.execute(event.context) == ProcessCode::SUCCESS);
  BOOST_CHECK(converter.collections() ==
              std::vector<std::string>({"CaloHits", "CaloLinks"}));
  const auto& hits =
      event.readPodio<edm4hep::CalorimeterHitCollection>(config.outputHits);
  const auto& links = event.readPodio<edm4hep::CaloHitSimCaloHitLinkCollection>(
      config.outputLinks);
  BOOST_REQUIRE_EQUAL(hits.size(), 2);
  BOOST_REQUIRE_EQUAL(links.size(), 3);
  BOOST_CHECK_EQUAL(hits[0].getCellID(), 42);
  BOOST_CHECK_EQUAL(hits[1].getCellID(), CalorimeterFixture::highCellId);
  BOOST_CHECK_CLOSE(hits[1].getEnergy(), 0.08, 1e-4);
  BOOST_CHECK_CLOSE(hits[1].getTime(), 4, 1e-4);
  BOOST_CHECK_EQUAL(hits[1].getPosition().x, 1500);
  double weightSum = 0;
  for (const auto& link : links) {
    BOOST_CHECK(link.getFrom().isAvailable());
    BOOST_CHECK(link.getTo().isAvailable());
    BOOST_CHECK_EQUAL(link.getFrom().getCellID(), link.getTo().getCellID());
    if (link.getFrom() == hits[1]) {
      const auto expected = link.getTo().getObjectID().index == 0 ? 0.75 : 0.25;
      BOOST_CHECK_CLOSE(link.getWeight(), expected, 1e-4);
      weightSum += link.getWeight();
    } else {
      BOOST_CHECK_EQUAL(link.getWeight(), 1);
    }
  }
  BOOST_CHECK_CLOSE(weightSum, 1, 1e-4);
  // Writing a relation does not consume or recalibrate native inputs.
  const auto& native = event.read<CalorimeterHitContainer>(config.inputHits);
  BOOST_CHECK_CLOSE(native[1].energy / 1_GeV, 0.08, 1e-4);
  BOOST_CHECK(native[1].sourceIndices == std::vector<std::size_t>({0, 1, 3}));
}

BOOST_AUTO_TEST_CASE(EmptyCollectionsAndZeroEnergyHits) {
  for (bool zeroHit : {false, true}) {
    podio::Frame frame;
    if (zeroHit) {
      frame = oneHit(0, {});
    } else {
      frame.put(edm4hep::SimCalorimeterHitCollection{}, "SimCaloHits");
    }
    Event event(std::move(frame));
    event.convert();
    BOOST_CHECK(
        event.read<SimCalorimeterHitContainer>("calo_deposits").empty());
    BOOST_CHECK(EDM4hepCalorimeterOutputConverter{outputConfig()}.execute(
                    event.context) == ProcessCode::SUCCESS);
    BOOST_CHECK(
        event.readPodio<edm4hep::CalorimeterHitCollection>("CaloHits").empty());
    BOOST_CHECK(
        event.readPodio<edm4hep::CaloHitSimCaloHitLinkCollection>("CaloLinks")
            .empty());
  }
}

BOOST_AUTO_TEST_CASE(RejectMissingOrWrongCollectionAndUnknownGeometry) {
  Event event;
  auto config = inputConfig();
  config.inputSimHits = "Missing";
  BOOST_CHECK_THROW(
      EDM4hepCalorimeterInputConverter{config}.execute(event.context),
      std::runtime_error);
  config.inputSimHits = "MCParticles";
  BOOST_CHECK_THROW(
      EDM4hepCalorimeterInputConverter{config}.execute(event.context),
      std::runtime_error);
  config = inputConfig();
  config.cellCentre = [](std::uint64_t) -> Acts::Vector3 {
    throw std::out_of_range("Unknown cell");
  };
  BOOST_CHECK_THROW(
      EDM4hepCalorimeterInputConverter{config}.execute(event.context),
      std::out_of_range);
  config.cellCentre = [](std::uint64_t) {
    return Acts::Vector3{std::numeric_limits<double>::quiet_NaN(), 0, 0};
  };
  BOOST_CHECK_THROW(
      EDM4hepCalorimeterInputConverter{config}.execute(event.context),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(RejectMalformedHitsContributionsAndMissingTimes) {
  const auto nan = std::numeric_limits<float>::quiet_NaN();
  const auto infinity = std::numeric_limits<float>::infinity();
  for (const auto energy : {-1.f, nan, infinity, 1.f}) {
    Event event(oneHit(energy, {}));
    BOOST_CHECK_THROW(
        EDM4hepCalorimeterInputConverter{inputConfig()}.execute(event.context),
        std::invalid_argument);
  }
  for (const auto values :
       {std::pair{-1.f, 1.f}, std::pair{nan, 1.f}, std::pair{infinity, 1.f},
        std::pair{1.f, nan}, std::pair{1.f, infinity}, std::pair{0.5f, 1.f}}) {
    Event event(oneHit(1, {values}));
    BOOST_CHECK_THROW(
        EDM4hepCalorimeterInputConverter{inputConfig()}.execute(event.context),
        std::invalid_argument);
  }
  Event repeated(oneHit(2, {{1, 1}}, true));
  BOOST_CHECK_THROW(
      EDM4hepCalorimeterInputConverter{inputConfig()}.execute(repeated.context),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(RejectMalformedCalibratedHitsAndProvenance) {
  const std::vector<std::function<void(CalorimeterHitContainer&)>> changes{
      [](auto& hits) { hits[0].sourceIndices.clear(); },
      [](auto& hits) { hits[0].sourceIndices = {999}; },
      [](auto& hits) { hits[0].sourceIndices = {4, 4}; },
      [](auto& hits) { hits[0].sourceIndices = {0}; },
      [](auto& hits) { hits[1].cellId = hits[0].cellId; },
      [](auto& hits) { hits[0].energy = -1; },
      [](auto& hits) {
        hits[0].time = std::numeric_limits<double>::infinity();
      },
      [](auto& hits) {
        hits[0].position.x() = std::numeric_limits<double>::quiet_NaN();
      }};
  for (const auto& change : changes) {
    Event event;
    event.convert();
    auto hits = event.read<CalorimeterHitContainer>("calo_hits");
    change(hits);
    event.write("changed_hits", std::move(hits));
    auto config = outputConfig();
    config.inputHits = "changed_hits";
    BOOST_CHECK_THROW(
        EDM4hepCalorimeterOutputConverter{config}.execute(event.context),
        std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(RejectSourceMappingMismatch) {
  for (int change = 0; change < 3; ++change) {
    Event event;
    event.convert();
    auto sources =
        event.read<EDM4hepCalorimeterSourceContainer>("calo_sources");
    if (change == 0) {
      sources.pop_back();
    } else if (change == 1) {
      sources[0].hit = sources[3].hit;
    } else {
      sources[0].contribution = sources[4].contribution;
    }
    event.write("changed_sources", std::move(sources));
    auto config = outputConfig();
    config.inputSources = "changed_sources";
    BOOST_CHECK_THROW(
        EDM4hepCalorimeterOutputConverter{config}.execute(event.context),
        std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(RejectFloatOverflowAndPositiveEnergyUnderflow) {
  for (int change = 0; change < 4; ++change) {
    Event event;
    event.convert();
    auto hits = event.read<CalorimeterHitContainer>("calo_hits");
    if (change == 0) {
      hits[0].energy = 2.0 * std::numeric_limits<float>::max() * 1_GeV;
    } else if (change == 1) {
      hits[0].energy = std::numeric_limits<double>::min() * 1_GeV;
    } else if (change == 2) {
      hits[0].position.x() = std::numeric_limits<double>::max();
    } else {
      hits[0].time = std::numeric_limits<double>::max();
    }
    event.write("changed_hits", std::move(hits));
    auto config = outputConfig();
    config.inputHits = "changed_hits";
    BOOST_CHECK_THROW(
        EDM4hepCalorimeterOutputConverter{config}.execute(event.context),
        std::overflow_error);
  }
}

BOOST_AUTO_TEST_CASE(RejectInvalidConverterConfiguration) {
  auto input = inputConfig();
  input.cellCentre = {};
  BOOST_CHECK_THROW(EDM4hepCalorimeterInputConverter{input},
                    std::invalid_argument);
  input = inputConfig();
  input.outputDeposits = input.outputSources;
  BOOST_CHECK_THROW(EDM4hepCalorimeterInputConverter{input},
                    std::invalid_argument);
  input = inputConfig();
  input.outputDeposits = input.inputFrame;
  BOOST_CHECK_THROW(EDM4hepCalorimeterInputConverter{input},
                    std::invalid_argument);
  input = inputConfig();
  input.inputSimHits.clear();
  BOOST_CHECK_THROW(EDM4hepCalorimeterInputConverter{input},
                    std::invalid_argument);
  auto output = outputConfig();
  output.outputHits.clear();
  BOOST_CHECK_THROW(EDM4hepCalorimeterOutputConverter{output},
                    std::invalid_argument);
  output = outputConfig();
  output.outputLinks = output.inputHits;
  BOOST_CHECK_THROW(EDM4hepCalorimeterOutputConverter{output},
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
