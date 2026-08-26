// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Mille/MillePedeSteering.hpp"

#include <filesystem>

using namespace ActsPlugins::ActsToMille;

BOOST_AUTO_TEST_SUITE(MillePedeSteeringTests)

/// catch a missing steering file
BOOST_AUTO_TEST_CASE(invalidSteeringDest) {
  MillePedeSteering::Config steerCfg;
  MillePedeSteering steer;
  auto res =
      steer.generateSteeringFile("/invalid/location/wontWork.txt", steerCfg);
  BOOST_CHECK(res.empty());
}

/// write a valid file
BOOST_AUTO_TEST_CASE(validFileCheck) {
  MillePedeSteering::Config steerCfg;
  MillePedeSteering steer;
  const std::filesystem::path testSteer = "testSteer.txt";
  auto res = steer.generateSteeringFile(testSteer, steerCfg);
  BOOST_CHECK(res == testSteer);
  BOOST_CHECK(std::filesystem::exists(testSteer) &&
              std::filesystem::is_regular_file(testSteer));
}

BOOST_AUTO_TEST_SUITE_END()
