// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <iomanip>

#include <boost/type_erasure/any_cast.hpp>

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include "ATLASCuts.hpp"
#include "SpacePoint.hpp"

#include "Acts/Utilities/Platforms/PlatformDef.h"
#include "Acts/Utilities/Platforms/CUDA/CudaUtils.cu"

std::vector<const SpacePoint*> readFile(std::string filename) {
  std::string line;
  int layer;
  std::vector<const SpacePoint*> readSP;

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
        SpacePoint* sp =
            new SpacePoint{x, y, z, r, layer, varianceR, varianceZ};
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
  auto start_pre = std::chrono::system_clock::now();
    
  std::string file{"sp.txt"};
  bool help(false);
  bool quiet(false);
  int  nGroupToIterate = 500;
  int  skip = 0;
  int  deviceID = 0;
    
  int opt;
  while ((opt = getopt(argc, argv, "hf:n:s:d:q")) != -1) {
    switch (opt) {
      case 'f':
        file = optarg;
        break;
      case 'n':
        nGroupToIterate = atoi(optarg);
        break;
      case 's':
        skip = atoi(optarg);
        break;		
      case 'd':
        deviceID = atoi(optarg);
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
              << "      -f FILE  : read spacepoints from FILE. Default is \""
              << file << "\"" << std::endl;
          std::cout
              << "      -n NUM   : Number of groups to iterate in seed finding. Default is "
              << nGroupToIterate << std::endl;
	  std::cout
              << "      -s SKIP  : Number of groups to skip in seed finding. Default is "
              << skip << std::endl;
          std::cout
              << "      -d DEVID : NVIDIA GPU device ID. Default is "
              << deviceID << std::endl;	  	  
          std::cout << "      -q : don't print out all found seeds"
                    << std::endl;
        }

        exit(EXIT_FAILURE);
    }
  }

  std::string devName;
  SetDevice(deviceID,devName);
  
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

  std::cout << "Number of groups to iterate: " << nGroupToIterate << std::endl;
  
  /// For CPU seed finder
  Acts::SeedfinderConfig<SpacePoint> config;
  // silicon detector max
  config.rMax = 160.;
  config.deltaRMin = 5.;
  config.deltaRMax = 160.;
  config.collisionRegionMin = -250.;
  config.collisionRegionMax = 250.;
  config.zMin = -2800.;
  config.zMax = 2800.;
  config.maxSeedsPerSpM = 5;
  // 2.7 eta
  config.cotThetaMax = 7.40627;
  config.sigmaScattering = 1.00000;

  config.minPt = 500.;
  config.bFieldInZ = 0.00199724;

  config.beamPos = {-.5, -.5};
  config.impactMax = 10.;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(sfconf, &atlasCuts));  
  Acts::Seedfinder<SpacePoint, Acts::CPU> seedfinder_cpu(config);
  Acts::Seedfinder<SpacePoint, Acts::CUDA> seedfinder_cuda(config);

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SpacePoint& sp, float, float, float) -> Acts::Vector2D {
    return {sp.varianceR, sp.varianceZ};
  };
  
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
  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConf);
  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(), spVec.end(), ct,
                                                 bottomBinFinder, topBinFinder,
                                                 std::move(grid), config);

  auto end_pre = std::chrono::system_clock::now();
  
  std::chrono::duration<double> elapsec_pre = end_pre - start_pre;
  std::cout << "Preprocess Time: " << elapsec_pre.count() << std::endl;
  
  int group_count;
  
  // CPU
  group_count=0;
  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector_cpu;
  auto groupIt = spGroup.begin();

  for (int i_s=0; i_s<skip; i_s++) ++groupIt;
  for (; !(groupIt == spGroup.end()); ++groupIt) {
    seedVector_cpu.push_back(seedfinder_cpu.createSeedsForGroup(
        groupIt.bottom(), groupIt.middle(), groupIt.top()));
    group_count++;
    if (group_count >= nGroupToIterate) break;
  }
  
  auto timeMetric_cpu = seedfinder_cpu.getTimeMetric();
  
  // CUDA
  //cudaProfilerStart();
  group_count=0;
  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector_cuda;
  groupIt = spGroup.begin();

  for (int i_s=0; i_s<skip; i_s++) ++groupIt;
  for (; !(groupIt == spGroup.end()); ++groupIt) {
    seedVector_cuda.push_back(seedfinder_cuda.createSeedsForGroup(
        groupIt.bottom(), groupIt.middle(), groupIt.top()));
    group_count++;
    if (group_count >= nGroupToIterate) break;
  } 

  auto timeMetric_cuda = seedfinder_cuda.getTimeMetric();


  std::cout << "------------Time Metric-------------" << std::endl;
  std::cout << "               CPU        CUDA  " << std::endl;
  std::cout << "DS time:   "
	    << std::setw(10) << std::get<0>(timeMetric_cpu)  << "  "
	    << std::setw(10) << std::get<0>(timeMetric_cuda) << std::endl;
  std::cout << "TC time:   "
	    << std::setw(10) << std::setw(10) << std::get<1>(timeMetric_cpu) << "  "
	    << std::setw(10) << std::get<1>(timeMetric_cuda) << std::endl;
  std::cout << "TS time:   "
	    << std::setw(10) << std::get<2>(timeMetric_cpu) << "  "
	    << std::setw(10) << std::get<2>(timeMetric_cuda) << std::endl;
  std::cout << "Wall time: "
	    << std::setw(10) << std::get<3>(timeMetric_cpu) << "  "
	    << std::setw(10) << std::get<3>(timeMetric_cuda) << std::endl;
  
  //cudaProfilerStop();
 
  int nSeed_cpu = 0;
  for (auto& outVec : seedVector_cpu) {
    nSeed_cpu += outVec.size();
  }

  int nSeed_cuda = 0;
  for (auto& outVec : seedVector_cuda) {
    nSeed_cuda += outVec.size();
  }
  
  std::cout << "Number of Seeds (CPU | CUDA): " << nSeed_cpu << " | " << nSeed_cuda << std::endl;

  int nMatch = 0;
  
  for (size_t i =0; i < seedVector_cpu.size(); i++){    
    auto regionVec_cpu  = seedVector_cpu[i];
    auto regionVec_cuda = seedVector_cuda[i];

    std::vector< std::vector< SpacePoint > > seeds_cpu;
    std::vector< std::vector< SpacePoint > > seeds_cuda;
              
    //for (size_t i_cpu = 0; i_cpu < regionVec_cpu.size(); i_cpu++) {
    for (auto sd: regionVec_cpu){
      std::vector< SpacePoint > seed_cpu;
      seed_cpu.push_back(*(sd.sp()[0]));
      seed_cpu.push_back(*(sd.sp()[1]));
      seed_cpu.push_back(*(sd.sp()[2]));

      seeds_cpu.push_back(seed_cpu);
    }

    for (auto sd: regionVec_cuda) {
      std::vector< SpacePoint > seed_cuda;           
      seed_cuda.push_back(*(sd.sp()[0]));
      seed_cuda.push_back(*(sd.sp()[1]));
      seed_cuda.push_back(*(sd.sp()[2]));

      seeds_cuda.push_back(seed_cuda);
    }
    
    for (auto seed: seeds_cpu){
      for (auto other: seeds_cuda){
	if (seed[0] == other[0] &&
	    seed[1] == other[1] &&
	    seed[2] == other[2]){
	  nMatch++;
	  break;
	}
      }
    }
  }

  std::cout << nMatch << " seeds are matched" << std::endl;
    
  if (!quiet) {
    std::cout << "CPU Seed result:" << std::endl;
    
    for (auto& regionVec : seedVector_cpu) {
      for (size_t i = 0; i < regionVec.size(); i++) {
        const Acts::Seed<SpacePoint>* seed = &regionVec[i];
        const SpacePoint* sp = seed->sp()[0];
        std::cout << " (" << sp->x() << ", " << sp->y() << ", " << sp->z()
                  << ") ";
        sp = seed->sp()[1];
        std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        sp = seed->sp()[2];
        std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        std::cout << std::endl;
      }
    }
    
    std::cout << std::endl;
    std::cout << "CUDA Seed result:" << std::endl;

    for (auto& regionVec : seedVector_cuda) {
      for (size_t i = 0; i < regionVec.size(); i++) {
        const Acts::Seed<SpacePoint>* seed = &regionVec[i];
        const SpacePoint* sp = seed->sp()[0];
        std::cout << " (" << sp->x() << ", " << sp->y() << ", " << sp->z()
                  << ") ";
        sp = seed->sp()[1];
        std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        sp = seed->sp()[2];
        std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        std::cout << std::endl;
      }
    }
  }

 
  return 0;
}


