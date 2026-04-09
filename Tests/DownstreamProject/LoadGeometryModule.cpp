// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryModuleLoader.hpp"

#ifdef ACTS_HAVE_DD4HEP
#include "ActsPlugins/DD4hep/GeometryModuleLoader.hpp"

#include <DD4hep/Detector.h>
#endif

#include <cstdlib>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  std::string modulePath;
#ifdef ACTS_HAVE_DD4HEP
  std::string dd4hepXml;
#endif

  po::options_description desc("Options");
  // clang-format off
  desc.add_options()
    ("help,h", "Show help")
    ("module", po::value<std::string>(&modulePath)->required(),
     "Path to geometry module shared library")
#ifdef ACTS_HAVE_DD4HEP
    ("dd4hep", po::value<std::string>(&dd4hepXml),
     "Path to DD4hep compact XML file; activates DD4hep detector loading")
#endif
    ;
  // clang-format on

  po::positional_options_description pos;
  pos.add("module", 1);

  po::variables_map vm;
  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pos).run(),
        vm);
    if (vm.count("help") != 0u) {
      std::cout << desc << std::endl;
      return EXIT_SUCCESS;
    }
    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << "Error: " << e.what() << "\n" << desc << std::endl;
    return EXIT_FAILURE;
  }

  try {
    auto logger =
        Acts::getDefaultLogger("LoadGeometryModule", Acts::Logging::VERBOSE);

    std::shared_ptr<Acts::TrackingGeometry> trackingGeometry;

#ifdef ACTS_HAVE_DD4HEP
    if (vm.count("dd4hep") != 0u) {
      auto detector = dd4hep::Detector::make_unique("LoadGeometryModule");
      detector->fromCompact(dd4hepXml);
      detector->volumeManager();
      detector->apply("DD4hepVolumeManager", 0, nullptr);
      trackingGeometry =
          Acts::loadDD4hepGeometryModule(modulePath, *detector, *logger);
    } else
#endif
    {
      trackingGeometry = Acts::loadGeometryModule(modulePath, *logger);
    }

    std::cout << "Geometry module loaded successfully" << std::endl;
    if (!trackingGeometry) {
      std::cerr << "Geometry module returned an invalid handle" << std::endl;
      return EXIT_FAILURE;
    }
  } catch (const std::exception& e) {
    std::cerr << "Unexpected failure while loading geometry module: "
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
