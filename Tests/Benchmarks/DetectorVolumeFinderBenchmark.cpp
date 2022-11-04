// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Common.hpp"
#include "Acts/Experimental/CylindricalDetectorHelper.hpp"
#include "Acts/Experimental/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Json/DetectorJsonConverter.hpp"
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>
#include <random>

#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace Acts;

int main(int argc, char* argv[]) {
  unsigned int toys = 1;
  unsigned int lvl = Acts::Logging::INFO;
  ActsScalar r = 1150.;
  ActsScalar hz = 2950.;

  // Create a test context
  GeometryContext tgContext = GeometryContext();

  try {
    po::options_description desc("Allowed options");
    // clang-format off
  desc.add_options()
      ("help", "produce help message")
      ("toys",po::value<unsigned int>(&toys)->default_value(20000),"number of tracks to propagate")
      ("verbose",po::value<unsigned int>(&lvl)->default_value(Acts::Logging::INFO),"logging level");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help") != 0u) {
      std::cout << desc << std::endl;
      return 0;
    }
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  ACTS_LOCAL_LOGGER(
      getDefaultLogger("DetectorVolumeFinder", Acts::Logging::Level(lvl)));

  // Read the detector
  // Let's assign the senstive surfaces
  auto oddIn = std::ifstream("ODD.json");
  nlohmann::json jodd;
  oddIn >> jodd;

  auto detector = Acts::Experimental::detectorFromJson(jodd);

  // Choose a random mean between 1 and 6
  std::mt19937 e0(1234);
  std::mt19937 e1(1234);
  std::mt19937 e2(1234);
  std::uniform_real_distribution<ActsScalar> r_dist(0, r);
  std::uniform_real_distribution<ActsScalar> phi_dist(-M_PI, M_PI);
  std::uniform_real_distribution<ActsScalar> z_dist(-hz, hz);

  size_t volumesFound = 0;
  size_t volumesMissed = 0;

  auto runBenchmark = [&](std::mt19937& egen,
                          const Experimental::Detector* detector,
                          const std::string& name) -> void {
    volumesFound = 0;
    volumesMissed = 0;
    // run the benchmark
    const auto benchmark = Acts::Test::microBenchmark(
        [&] {
          auto r = r_dist(egen);
          auto phi = phi_dist(egen);
          auto z = z_dist(egen);
          Vector3 position(r * std::cos(phi), r * std::sin(phi), z);
          if (detector != nullptr) {
            auto volume = detector->findVolume(tgContext, position);
            if (volume != nullptr) {
              ++volumesFound;
            } else {
              ++volumesMissed;
            }
          }
        },
        1, toys);

    ACTS_INFO("Test benchmark: " << name);
    ACTS_INFO("                 - execution stats: " << benchmark);
    if (detector != nullptr) {
      ACTS_INFO("                 - hit/missed     : " << volumesFound << "/"
                                                       << volumesMissed);
    }
  };

  runBenchmark(e0, nullptr, "empty");
  runBenchmark(e1, detector.get(), "trial&error");
  Acts::Experimental::attachGridVolumeFinder(tgContext, detector,
                                             Acts::Logging::INFO);
  runBenchmark(e1, detector.get(), "grid");
  // 
  return 0;
}
