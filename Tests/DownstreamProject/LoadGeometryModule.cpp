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

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <module-path>" << std::endl;
    return EXIT_FAILURE;
  }

  try {
    auto logger =
        Acts::getDefaultLogger("LoadGeometryModule", Acts::Logging::VERBOSE);
    auto trackingGeometry = Acts::loadGeometryModule(argv[1], *logger);
    std::cout << "Geometry module loaded successfully" << std::endl;
    if (!trackingGeometry) {
      std::cerr << "Geometry module returned an invalid handle" << std::endl;
      return EXIT_FAILURE;
    }
  } catch (const std::exception& e) {
    std::cerr << "Unexpected failure while loading good module: " << e.what()
              << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
