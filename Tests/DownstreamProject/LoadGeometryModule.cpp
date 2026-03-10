// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryModuleLoader.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string_view>

int main() {
  auto request = Acts::makeGeometryModuleRequest();

  try {
    auto handle = Acts::loadGeometryModule(ACTS_TEST_GEOMETRY_MODULE_PATH,
                                           request);
    if (!handle) {
      std::cerr << "Geometry module returned an invalid handle" << std::endl;
      return EXIT_FAILURE;
    }
  } catch (const std::exception& e) {
    std::cerr << "Unexpected failure while loading good module: " << e.what()
              << std::endl;
    return EXIT_FAILURE;
  }

  try {
    auto mismatchedRequest = request;
    mismatchedRequest.abi_tag = "intentionally-mismatched-abi-tag";
    auto mismatchHandle = Acts::loadGeometryModule(ACTS_TEST_GEOMETRY_MODULE_PATH,
                                                   mismatchedRequest);
    (void)mismatchHandle;
    std::cerr << "Expected ABI mismatch for mismatched request" << std::endl;
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::string_view message(e.what());
    if (message.find("ABI mismatch") == std::string_view::npos) {
      std::cerr << "Mismatch module failed for unexpected reason: " << e.what()
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
