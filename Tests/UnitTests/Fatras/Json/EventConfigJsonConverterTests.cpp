// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/PdgParticle.hpp"
#include "ActsFatras/Json/EventConfigJsonConverter.hpp"
#include "ActsFatras/Synthetic/EventConfig.hpp"

#include <filesystem>
#include <fstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

using namespace ActsFatras::Synthetic;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(SyntheticEventConfigJsonSuite)

/// A configuration survives a round trip, field for field.
BOOST_AUTO_TEST_CASE(EventConfigRoundTrip) {
  EventConfig original;
  // moved off every default, so that a field written or read as the wrong one
  // cannot pass unnoticed
  original.generation.pileup = 37;
  original.generation.chargedPerUnitRapidity = 4.25f;
  original.generation.minPt = 0.031f;
  original.generation.rapidityEdgeWidth = 0.17f;
  original.simulation.propagation.maxTurns = 0.75f;
  original.simulation.material.scale = 1.3f;
  original.simulation.material.energyLossModel = EnergyLossModel::Mean;
  original.simulation.material.multipleScattering = false;
  original.simulation.measurement.positionSmearing = 0.011f;
  original.simulation.secondaries.electronRate = 3.75f;
  original.simulation.secondaries.maxGenerations = 3;
  original.simulation.secondaries.sampling.kt = 0.271f;
  original.particlePdg = Acts::PdgParticle::eMuon;
  original.bFieldZ = 1.997f;
  original.seed = 4242;

  const nlohmann::json written = original;
  const EventConfig read = written.get<EventConfig>();
  BOOST_CHECK(nlohmann::json(read) == written);

  BOOST_CHECK_EQUAL(read.generation.pileup, 37u);
  BOOST_CHECK_EQUAL(read.generation.chargedPerUnitRapidity, 4.25f);
  BOOST_CHECK_EQUAL(read.simulation.propagation.maxTurns, 0.75f);
  BOOST_CHECK_EQUAL(read.simulation.material.scale, 1.3f);
  BOOST_CHECK(read.simulation.material.energyLossModel ==
              EnergyLossModel::Mean);
  BOOST_CHECK(!read.simulation.material.multipleScattering);
  BOOST_CHECK_EQUAL(read.simulation.secondaries.maxGenerations, 3u);
  BOOST_CHECK_EQUAL(read.simulation.secondaries.sampling.kt, 0.271f);
  BOOST_CHECK(read.particlePdg == Acts::PdgParticle::eMuon);
  BOOST_CHECK_EQUAL(read.seed, 4242u);
}

/// Through a file, header and all.
BOOST_AUTO_TEST_CASE(EventConfigRoundTripsThroughAFile) {
  EventConfig original;
  original.generation.pileup = 5;
  original.seed = 7;

  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "acts-synthetic-config.json";
  writeEventConfig(path, original);
  const EventConfig read = readEventConfig(path);
  BOOST_CHECK(nlohmann::json(read) == nlohmann::json(original));
  std::filesystem::remove(path);

  const std::filesystem::path other =
      std::filesystem::temp_directory_path() / "acts-synthetic-no-header.json";
  {
    std::ofstream file(other);
    file << nlohmann::json(original).dump(2) << std::endl;
  }
  BOOST_CHECK_THROW(readEventConfig(other), std::runtime_error);
  std::filesystem::remove(other);
}

/// A configuration is fitted as a whole, so a field left out of a file is an
/// error rather than a default: taking the default would retune every other
/// number that was fitted alongside it.
BOOST_AUTO_TEST_CASE(EventConfigNeedsEveryField) {
  const nlohmann::json complete = EventConfig{};
  BOOST_CHECK_NO_THROW(complete.get<EventConfig>());

  for (const char* section : {"generation", "simulation"}) {
    nlohmann::json missing = complete;
    missing.erase(section);
    BOOST_CHECK_THROW(missing.get<EventConfig>(), nlohmann::json::out_of_range);
  }

  nlohmann::json missingLeaf = complete;
  missingLeaf.at("simulation").at("secondaries").erase("electronRate");
  BOOST_CHECK_THROW(missingLeaf.get<EventConfig>(),
                    nlohmann::json::out_of_range);

  nlohmann::json missingDeepLeaf = complete;
  missingDeepLeaf.at("simulation").at("secondaries").at("sampling").erase("kt");
  BOOST_CHECK_THROW(missingDeepLeaf.get<EventConfig>(),
                    nlohmann::json::out_of_range);
}

/// An enum is written as what it means, not as the integer it happens to be.
BOOST_AUTO_TEST_CASE(EnergyLossModelIsNamed) {
  const nlohmann::json config = EventConfig{};
  BOOST_CHECK_EQUAL(
      config.at("simulation").at("material").at("energyLossModel"), "mode");

  EventConfig mean;
  mean.simulation.material.energyLossModel = EnergyLossModel::Mean;
  BOOST_CHECK_EQUAL(nlohmann::json(mean)
                        .at("simulation")
                        .at("material")
                        .at("energyLossModel"),
                    "mean");

  nlohmann::json wrong = config;
  wrong["simulation"]["material"]["energyLossModel"] = "most-probable";
  BOOST_CHECK_THROW(wrong.get<EventConfig>(), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
