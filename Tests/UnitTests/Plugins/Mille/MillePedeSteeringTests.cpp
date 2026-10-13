// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Mille/MillePedeError.hpp"
#include "ActsPlugins/Mille/MillePedeSteering.hpp"

#include <filesystem>

using namespace ActsPlugins;

BOOST_AUTO_TEST_SUITE(MillePedeSteeringTests)

/// catch a missing steering file
BOOST_AUTO_TEST_CASE(invalidSteeringDest) {
  auto logger =
      Acts::getDefaultLogger("invalidSteeringDest", Acts::Logging::INFO);
  // This test deliberately emits an error
  // message. To prevent the ACTS CI from
  // interpreting this as a failure,
  // we temporarily disable the failure threshold
  Acts::Logging::ScopedFailureThreshold threshold{Acts::Logging::Level::MAX};
  MillePedeSteeringConfig steerCfg;
  auto res = ActsPlugins::generateMillePedeSteeringFile(
      "/invalid/location/wontWork.txt", steerCfg, *logger);
  BOOST_CHECK(!res.ok());
  BOOST_CHECK_EQUAL(res.error(), MillePedeError::UnableToWriteSteering);
}

/// write a valid file
BOOST_AUTO_TEST_CASE(validFileCheck) {
  auto logger = Acts::getDefaultLogger("validFileCheck", Acts::Logging::INFO);
  MillePedeSteeringConfig steerCfg{.inputFiles = {"DummyBinary.dat"}};
  const std::filesystem::path testSteer = "testSteer.txt";
  auto res =
      ActsPlugins::generateMillePedeSteeringFile(testSteer, steerCfg, *logger);
  BOOST_CHECK(res.ok());
  BOOST_CHECK(res.value() == testSteer);
  BOOST_CHECK(std::filesystem::exists(testSteer) &&
              std::filesystem::is_regular_file(testSteer));
}

BOOST_AUTO_TEST_SUITE_END()
