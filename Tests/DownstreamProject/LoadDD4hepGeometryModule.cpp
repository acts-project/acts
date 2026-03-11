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
#include <memory>
#include <stdexcept>

#include <DD4hep/Detector.h>

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <dd4hep-xml> <module-path>"
              << std::endl;
    return EXIT_FAILURE;
  }

  const char* xmlPath = argv[1];
  const char* modulePath = argv[2];

  try {
    auto logger = Acts::getDefaultLogger("LoadDD4hepGeometryModule",
                                         Acts::Logging::VERBOSE);

    auto detector = dd4hep::Detector::make_unique("LoadDD4hepGeometryModule");
    detector->fromCompact(xmlPath);
    detector->volumeManager();
    detector->apply("DD4hepVolumeManager", 0, nullptr);

    auto trackingGeometry =
        Acts::loadGeometryModule(modulePath, &*detector, *logger);

    std::cout << "DD4hep geometry module loaded successfully" << std::endl;
    if (!trackingGeometry) {
      std::cerr << "Geometry module returned an invalid handle" << std::endl;
      return EXIT_FAILURE;
    }
  } catch (const std::exception& e) {
    std::cerr << "Unexpected failure while loading DD4hep geometry module: "
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
