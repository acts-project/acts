// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePoint.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

#include <boost/type_erasure/any_cast.hpp>

#include "ATLASCuts.hpp"

using namespace Acts::UnitLiterals;

class LayerLink : public Acts::SourceLink {
 public:
  LayerLink(size_t layer) : SourceLink(layer) {}
};

std::vector<Acts::SpacePoint*> readFile(std::string filename) {
  std::string line;
  size_t layer;
  std::vector<Acts::SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()) {
    while (!spFile.eof()) {
      std::getline(spFile, line);
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      float x, y, z, r, varianceR, varianceZ;
      if (linetype == "lxyz") {
        ss >> layer >> x >> y >> z >> varianceR >> varianceZ;
        r = std::sqrt(x * x + y * y);
        float f22 = varianceR;
        float wid = varianceZ;
        float cov = wid * wid * .08333;
        if (cov < f22)
          cov = f22;
        if (std::abs(z) > 450.) {
          varianceZ = 9. * cov;
          varianceR = .06;
        } else {
          varianceR = 9. * cov;
          varianceZ = .06;
        }

        LayerLink layerLink{layer};
        Acts::SourceLink* sourceLink =
            dynamic_cast<Acts::SourceLink*>(&layerLink);

        auto sp = new Acts::SpacePoint{Acts::Vector3(x, y, z),
                                       Acts::Vector2(-.5_mm, -.5_mm),
                                       Acts::Vector2(varianceR, varianceZ),
                                       Acts::Vector2(0., 0.),
                                       5.,
                                       {sourceLink}};
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

  int opt;
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
  std::vector<Acts::SpacePoint*> spVec = readFile(file);
  auto end_read = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_read = end_read - start_read;

  std::cout << "read " << spVec.size() << " SP from file " << file << " in "
            << elapsed_read.count() << "s" << std::endl;

  Acts::SeedfinderConfig config;
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
  config.bFieldInZ = 1.99724_T;

  config.beamPos = {-.5_mm, -.5_mm};
  config.impactMax = 10._mm;

  config.useVariableMiddleSPRange = false;
  Acts::Extent rRangeSPExtent;

  int numPhiNeighbors = 1;

  std::vector<std::pair<int, int>> zBinNeighborsTop;
  std::vector<std::pair<int, int>> zBinNeighborsBottom;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder>(
      Acts::BinFinder(zBinNeighborsBottom, numPhiNeighbors));
  auto topBinFinder = std::make_shared<Acts::BinFinder>(
      Acts::BinFinder(zBinNeighborsTop, numPhiNeighbors));
  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts atlasCuts = Acts::ATLASCuts();
  config.seedFilter =
      std::make_unique<Acts::SeedFilter>(Acts::SeedFilter(sfconf, &atlasCuts));
  Acts::Seedfinder seedfinder(config);

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;
  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid> grid =
      Acts::SpacePointGridCreator::createGrid(gridConf);
  auto spGroup =
      Acts::BinnedSPGroup(spVec.begin(), spVec.end(), bottomBinFinder,
                          topBinFinder, std::move(grid), config);

  std::vector<std::vector<Acts::InternalSeed>> seedVector;
  decltype(seedfinder)::State state;
  auto start = std::chrono::system_clock::now();
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();
  for (; !(groupIt == endOfGroups); ++groupIt) {
    auto& v = seedVector.emplace_back();
    seedfinder.createSeedsForGroup(state, std::back_inserter(v),
                                   groupIt.bottom(), groupIt.middle(),
                                   groupIt.top(), rRangeSPExtent);
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
      for (size_t i = 0; i < regionVec.size(); i++) {
        Acts::InternalSeed* seed = &regionVec[i];
        Acts::SpacePoint* sp = seed->sp[0];
        std::cout << " (" << sp->x() << ", " << sp->y() << ", " << sp->z()
                  << ") ";
        sp = seed->sp[1];
        std::cout << sp->getSourceLinks()[0]->geometryId() << " (" << sp->x()
                  << ", " << sp->y() << ", " << sp->z() << ") ";
        sp = seed->sp[2];
        std::cout << sp->getSourceLinks()[0]->geometryId() << " (" << sp->x()
                  << ", " << sp->y() << ", " << sp->z() << ") ";
        std::cout << std::endl;
      }
    }
  }
  return 0;
}
