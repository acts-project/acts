// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ATLASCuts.hpp"
#include "SpacePoint.hpp"

using namespace Acts::UnitLiterals;

std::vector<const SpacePoint*> readFile(const std::string& filename) {
  std::string line;
  int layer = 0;
  std::vector<const SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()) {
    while (!spFile.eof()) {
      std::getline(spFile, line);
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      if (linetype == "lxyz") {
        float x = 0, y = 0, z = 0, varianceR = 0, varianceZ = 0;
        std::optional<float> t, varianceT;
        ss >> layer >> x >> y >> z >> varianceR >> varianceZ;
        const float r = std::hypot(x, y);

        float cov = varianceZ * varianceZ * .08333;
        if (cov < varianceR) {
          cov = varianceR;
        }

        if (std::abs(z) > 450.) {
          varianceZ = 9. * cov;
          varianceR = .06;
        } else {
          varianceR = 9. * cov;
          varianceZ = .06;
        }

        SpacePoint* sp = new SpacePoint{
            x, y, z, r, layer, varianceR, varianceZ, t, varianceT};
        //     if(r < 200.){
        //       sp->setClusterList(1,0);
        //     }
        readSP.push_back(sp);
      }
    }
  }
  return readSP;
}

int main(int argc, char** argv) {
  std::string file{"sp.txt"};
  bool help(false);
  bool quiet(false);

  int opt = -1;
  while ((opt = getopt(argc, argv, "hf:q")) != -1) {
    switch (opt) {
      case 'f':
        file = optarg;
        break;
      case 'q':
        quiet = true;
        break;
      case 'h':
        help = true;
        [[fallthrough]];
      default: /* '?' */
        std::cerr << "Usage: " << argv[0] << " [-hq] [-f FILENAME]\n";
        if (help) {
          std::cout << "      -h : this help" << std::endl;
          std::cout
              << "      -f FILE : read spacepoints from FILE. Default is \""
              << file << "\"" << std::endl;
          std::cout << "      -q : don't print out all found seeds"
                    << std::endl;
        }

        exit(EXIT_FAILURE);
    }
  }

  std::ifstream f(file);
  if (!f.good()) {
    std::cerr << "input file \"" << file << "\" does not exist\n";
    exit(EXIT_FAILURE);
  }

  auto start_read = std::chrono::system_clock::now();
  std::vector<const SpacePoint*> spVec = readFile(file);
  auto end_read = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_read = end_read - start_read;

  std::cout << "read " << spVec.size() << " SP from file " << file << " in "
            << elapsed_read.count() << "s" << std::endl;

  Acts::SeedFinderConfig<SpacePoint> config;
  // silicon detector max
  config.rMax = 160._mm;
  config.deltaRMin = 5._mm;
  config.deltaRMax = 160._mm;
  config.deltaRMinTopSP = config.deltaRMin;
  config.deltaRMinBottomSP = config.deltaRMin;
  config.deltaRMaxTopSP = config.deltaRMax;
  config.deltaRMaxBottomSP = config.deltaRMax;
  config.collisionRegionMin = -250._mm;
  config.collisionRegionMax = 250._mm;
  config.zMin = -2800._mm;
  config.zMax = 2800._mm;
  config.maxSeedsPerSpM = 5;
  // 2.7 eta
  config.cotThetaMax = 7.40627;
  config.sigmaScattering = 1.00000;

  config.minPt = 500._MeV;

  config.impactMax = 10._mm;

  config.useVariableMiddleSPRange = false;

  Acts::SeedFinderOptions options;
  options.beamPos = {-.5_mm, -.5_mm};
  options.bFieldInZ = 2_T;

  int numPhiNeighbors = 1;

  // extent used to store r range for middle spacepoint
  Acts::Extent rRangeSPExtent;

  config.useVariableMiddleSPRange = false;
  const Acts::Range1D<float> rMiddleSPRange;

  std::vector<std::pair<int, int>> zBinNeighborsTop;
  std::vector<std::pair<int, int>> zBinNeighborsBottom;

  auto bottomBinFinder = std::make_unique<Acts::GridBinFinder<2ul>>(
      Acts::GridBinFinder<2ul>(numPhiNeighbors, zBinNeighborsBottom));
  auto topBinFinder = std::make_unique<Acts::GridBinFinder<2ul>>(
      Acts::GridBinFinder<2ul>(numPhiNeighbors, zBinNeighborsTop));
  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(sfconf, &atlasCuts));
  Acts::SeedFinder<SpacePoint, Acts::CylindricalSpacePointGrid<SpacePoint>>
      a;  // test creation of unconfigured finder
  a = Acts::SeedFinder<SpacePoint, Acts::CylindricalSpacePointGrid<SpacePoint>>(
      config);

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SpacePoint& sp, float, float, float) {
    Acts::Vector3 position(sp.x(), sp.y(), sp.z());
    Acts::Vector2 covariance(sp.varianceR, sp.varianceZ);
    return std::make_tuple(position, covariance, sp.t());
  };

  // setup spacepoint grid config
  Acts::CylindricalSpacePointGridConfig gridConf;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;
  // setup spacepoint grid options
  Acts::CylindricalSpacePointGridOptions gridOpts;
  gridOpts.bFieldInZ = options.bFieldInZ;
  // create grid with bin sizes according to the configured geometry

  Acts::CylindricalSpacePointGrid<SpacePoint> grid =
      Acts::CylindricalSpacePointGridCreator::createGrid<SpacePoint>(gridConf,
                                                                     gridOpts);
  Acts::CylindricalSpacePointGridCreator::fillGrid(
      config, options, grid, spVec.begin(), spVec.end(), ct, rRangeSPExtent);
  auto spGroup = Acts::CylindricalBinnedGroup<SpacePoint>(
      std::move(grid), *bottomBinFinder, *topBinFinder);

  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector;
  decltype(a)::SeedingState state;
  auto start = std::chrono::system_clock::now();
  for (auto [bottom, middle, top] : spGroup) {
    auto& v = seedVector.emplace_back();
    a.createSeedsForGroup(options, state, spGroup.grid(), std::back_inserter(v),
                          bottom, middle, top, rMiddleSPRange);
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "time to create seeds: " << elapsed_seconds.count() << std::endl;
  std::cout << "Number of regions: " << seedVector.size() << std::endl;
  int numSeeds = 0;
  for (auto& outVec : seedVector) {
    numSeeds += outVec.size();
  }
  std::cout << "Number of seeds generated: " << numSeeds << std::endl;
  if (!quiet) {
    for (auto& regionVec : seedVector) {
      for (std::size_t i = 0; i < regionVec.size(); i++) {
        const Acts::Seed<SpacePoint>* seed = &regionVec[i];
        const SpacePoint* sp = seed->sp()[0];
        std::cout << " (" << sp->x() << ", " << sp->y() << ", " << sp->z()
                  << ") ";
        sp = seed->sp()[1];
        std::cout << sp->layer << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        sp = seed->sp()[2];
        std::cout << sp->layer << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        std::cout << std::endl;
      }
    }
  }
  return 0;
}
