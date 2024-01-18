// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/Seed.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Plugins/Cuda/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/EventData/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <utility>

#include <boost/type_erasure/any_cast.hpp>
#include <cuda_profiler_api.h>

#include "ATLASCuts.hpp"
#include "SpacePoint.hpp"
#include "SpacePointContainer.hpp"

using namespace Acts::UnitLiterals;

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
        r = std::hypot(x, y);
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
  bool allgroup(false);
  bool do_cpu(true);
  int nGroupToIterate = 500;
  int skip = 0;
  int deviceID = 0;
  int nTrplPerSpBLimit = 100;
  int nAvgTrplPerSpBLimit = 2;

  int opt;
  while ((opt = getopt(argc, argv, "haf:n:s:d:l:m:qG")) != -1) {
    switch (opt) {
      case 'a':
        allgroup = true;
        break;
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
      case 'l':
        nAvgTrplPerSpBLimit = atoi(optarg);
        break;
      case 'm':
        nTrplPerSpBLimit = atoi(optarg);
        break;
      case 'q':
        quiet = true;
        break;
      case 'G':
        do_cpu = false;
        break;
      case 'h':
        help = true;
        [[fallthrough]];
      default: /* '?' */
        std::cerr << "Usage: " << argv[0] << " [-hq] [-f FILENAME]\n";
        if (help) {
          std::cout << "      -h : this help" << std::endl;
          std::cout << "      -a ALL   : analyze all groups. Default is \""
                    << allgroup << "\"" << std::endl;
          std::cout
              << "      -f FILE  : read spacepoints from FILE. Default is \""
              << file << "\"" << std::endl;
          std::cout << "      -n NUM   : Number of groups to iterate in seed "
                       "finding. Default is "
                    << nGroupToIterate << std::endl;
          std::cout << "      -s SKIP  : Number of groups to skip in seed "
                       "finding. Default is "
                    << skip << std::endl;
          std::cout << "      -d DEVID : NVIDIA GPU device ID. Default is "
                    << deviceID << std::endl;
          std::cout << "      -l : A limit on the average number of triplets "
                       "per bottom spacepoint: this is used for determining "
                       "matrix size for triplets per middle space point"
                    << nAvgTrplPerSpBLimit << std::endl;
          std::cout << "      -m : A limit on the number of triplets per "
                       "bottom spacepoint: users do not have to touch this for "
                       "# spacepoints < ~200k"
                    << nTrplPerSpBLimit << std::endl;
          std::cout << "      -q : don't print out all found seeds"
                    << std::endl;
          std::cout << "      -G : only run on GPU, not CPU" << std::endl;
        }

        exit(EXIT_FAILURE);
    }
  }

  std::string devName;
  ACTS_CUDA_ERROR_CHECK(cudaSetDevice(deviceID));

  std::ifstream f(file);
  if (!f.good()) {
    std::cerr << "input file \"" << file << "\" does not exist\n";
    exit(EXIT_FAILURE);
  }

  std::vector<const SpacePoint*> spVec = readFile(file);

  // Config
  Acts::SpacePointContainerConfig spConfig;
  // Options
  Acts::SpacePointContainerOptions spOptions;
  spOptions.beamPos = {-.5_mm, -.5_mm};

  // Prepare interface SpacePoint backend-ACTS
  ActsExamples::SpacePointContainer container(spVec);
  // Prepare Acts API
  Acts::SpacePointContainer<decltype(container), Acts::detail::RefHolder>
      spContainer(spConfig, spOptions, container);

  using value_type = typename decltype(spContainer)::ConstSpacePointProxyType;
  using seed_type = Acts::Seed<value_type>;

  std::cout << "read " << spContainer.size() << " SP from file " << file
            << std::endl;

  // Set seed finder configuration
  Acts::SeedFinderConfig<value_type> config;
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
  config = config.toInternalUnits();

  Acts::SeedFinderOptions options;
  options.bFieldInZ = 2_T;
  options.beamPos = {-.5_mm, -.5_mm};
  options = options.toInternalUnits().calculateDerivedQuantities(config);

  int numPhiNeighbors = 1;

  // extent used to store r range for middle spacepoint
  Acts::Extent rRangeSPExtent;

  const Acts::Range1D<float> rMiddleSPRange;

  std::vector<std::pair<int, int>> zBinNeighborsTop;
  std::vector<std::pair<int, int>> zBinNeighborsBottom;

  // cuda
  cudaDeviceProp prop;
  ACTS_CUDA_ERROR_CHECK(cudaGetDeviceProperties(&prop, deviceID));
  printf("\n GPU Device %d: \"%s\" with compute capability %d.%d\n\n", deviceID,
         prop.name, prop.major, prop.minor);
  config.maxBlockSize = prop.maxThreadsPerBlock;
  config.nTrplPerSpBLimit = nTrplPerSpBLimit;
  config.nAvgTrplPerSpBLimit = nAvgTrplPerSpBLimit;

  // binfinder
  auto bottomBinFinder = std::make_unique<Acts::GridBinFinder<2ul>>(
      Acts::GridBinFinder<2ul>(numPhiNeighbors, zBinNeighborsBottom));
  auto topBinFinder = std::make_unique<Acts::GridBinFinder<2ul>>(
      Acts::GridBinFinder<2ul>(numPhiNeighbors, zBinNeighborsTop));
  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts<value_type> atlasCuts = Acts::ATLASCuts<value_type>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<value_type>>(
      Acts::SeedFilter<value_type>(sfconf, &atlasCuts));
  Acts::SeedFinder<value_type> seedFinder_cpu(config);
  Acts::SeedFinder<value_type, Acts::Cuda> seedFinder_cuda(config, options);

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;
  // setup spacepoint grid options
  Acts::SpacePointGridOptions gridOpts;
  gridOpts.bFieldInZ = options.bFieldInZ;
  // create grid with bin sizes according to the configured geometry
  Acts::SpacePointGrid<value_type> grid =
      Acts::SpacePointGridCreator::createGrid<value_type>(gridConf, gridOpts);
  Acts::SpacePointGridCreator::fillGrid(config, options, grid, spVec.begin(),
                                        spVec.end(), ct, rRangeSPExtent);

  auto spGroup = Acts::BinnedSPGroup<value_type>(
      std::move(grid), *bottomBinFinder, *topBinFinder);

  auto end_pre = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsec_pre = end_pre - start_pre;
  double preprocessTime = elapsec_pre.count();
  std::cout << "Preprocess Time: " << preprocessTime << std::endl;

  //--------------------------------------------------------------------//
  //                        Begin Seed finding                          //
  //--------------------------------------------------------------------//

  auto start_cpu = std::chrono::system_clock::now();

  std::array<std::vector<std::size_t>, 2ul> navigation;
  navigation[0ul].resize(spGroup.grid().numLocalBins()[0ul]);
  navigation[1ul].resize(spGroup.grid().numLocalBins()[1ul]);
  std::iota(navigation[0ul].begin(), navigation[0ul].end(), 1ul);
  std::iota(navigation[1ul].begin(), navigation[1ul].end(), 1ul);

  std::array<std::size_t, 2ul> localPosition =
      spGroup.grid().localBinsFromGlobalBin(skip);

  int group_count;
  auto groupIt = Acts::BinnedSPGroupIterator<value_type>(spGroup, localPosition,
                                                         navigation);

  //----------- CPU ----------//
  group_count = 0;
  std::vector<std::vector<seed_type>> seedVector_cpu;

  if (do_cpu) {
    decltype(seedFinder_cpu)::SeedingState state;

    for (; groupIt != spGroup.end(); ++groupIt) {
      const auto [bottom, middle, top] = *groupIt;
      seedFinder_cpu.createSeedsForGroup(
          options, state, spGroup.grid(),
          std::back_inserter(seedVector_cpu.emplace_back()), bottom, middle,
          top, rMiddleSPRange);
      group_count++;
      if (allgroup == false) {
        if (group_count >= nGroupToIterate)
          break;
      }
    }
    // auto timeMetric_cpu = seedFinder_cpu.getTimeMetric();
    std::cout << "Analyzed " << group_count << " groups for CPU" << std::endl;
  }

  auto end_cpu = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsec_cpu = end_cpu - start_cpu;
  double cpuTime = elapsec_cpu.count();

  //----------- CUDA ----------//

  cudaProfilerStart();
  auto start_cuda = std::chrono::system_clock::now();

  group_count = 0;

  std::vector<std::vector<seed_type>> seedVector_cuda;
  groupIt = Acts::BinnedSPGroupIterator<value_type>(spGroup, localPosition,
                                                    navigation);

  Acts::SpacePointData spacePointData;
  spacePointData.resize(spVec.size());
>>>>>>> upstream/main

  for (; groupIt != spGroup.end(); ++groupIt) {
    const auto [bottom, middle, top] = *groupIt;
    seedVector_cuda.push_back(seedFinder_cuda.createSeedsForGroup(
        spGroup.grid(), bottom, middle, top));
    group_count++;
    if (allgroup == false) {
      if (group_count >= nGroupToIterate)
        break;
    }
  }

  auto end_cuda = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsec_cuda = end_cuda - start_cuda;
  double cudaTime = elapsec_cuda.count();

  cudaProfilerStop();
  std::cout << "Analyzed " << group_count << " groups for CUDA" << std::endl;

  std::cout << std::endl;
  std::cout << "----------------------- Time Metric -----------------------"
            << std::endl;
  std::cout << "                       " << (do_cpu ? "CPU" : "   ")
            << "          CUDA        " << (do_cpu ? "Speedup " : "")
            << std::endl;
  std::cout << "Seedfinding_Time  " << std::setw(11)
            << (do_cpu ? std::to_string(cpuTime) : "") << "  " << std::setw(11)
            << cudaTime << "  " << std::setw(11)
            << (do_cpu ? std::to_string(cpuTime / cudaTime) : "") << std::endl;
  double wallTime_cpu = cpuTime + preprocessTime;
  double wallTime_cuda = cudaTime + preprocessTime;
  std::cout << "Wall_time         " << std::setw(11)
            << (do_cpu ? std::to_string(wallTime_cpu) : "") << "  "
            << std::setw(11) << wallTime_cuda << "  " << std::setw(11)
            << (do_cpu ? std::to_string(wallTime_cpu / wallTime_cuda) : "")
            << std::endl;
  std::cout << "-----------------------------------------------------------"
            << std::endl;
  std::cout << std::endl;

  int nSeed_cpu = 0;
  for (auto& outVec : seedVector_cpu) {
    nSeed_cpu += outVec.size();
  }

  int nSeed_cuda = 0;
  for (auto& outVec : seedVector_cuda) {
    nSeed_cuda += outVec.size();
  }

  std::cout << "Number of Seeds (CPU | CUDA): " << nSeed_cpu << " | "
            << nSeed_cuda << std::endl;

  int nMatch = 0;

  for (std::size_t i = 0; i < seedVector_cpu.size(); i++) {
    auto regionVec_cpu = seedVector_cpu[i];
    auto regionVec_cuda = seedVector_cuda[i];

    std::vector<std::vector<value_type>> seeds_cpu;
    std::vector<std::vector<value_type>> seeds_cuda;

    // for (std::size_t i_cpu = 0; i_cpu < regionVec_cpu.size(); i_cpu++) {
    for (auto sd : regionVec_cpu) {
      std::vector<value_type> seed_cpu;
      seed_cpu.push_back(*(sd.sp()[0]));
      seed_cpu.push_back(*(sd.sp()[1]));
      seed_cpu.push_back(*(sd.sp()[2]));

      seeds_cpu.push_back(seed_cpu);
    }

    for (auto sd : regionVec_cuda) {
      std::vector<value_type> seed_cuda;
      seed_cuda.push_back(*(sd.sp()[0]));
      seed_cuda.push_back(*(sd.sp()[1]));
      seed_cuda.push_back(*(sd.sp()[2]));

      seeds_cuda.push_back(seed_cuda);
    }

    for (auto seed : seeds_cpu) {
      for (auto other : seeds_cuda) {
        if (*seed[0].sp() == *other[0].sp() and
            *seed[1].sp() == *other[1].sp() and
            *seed[2].sp() == *other[2].sp()) {
          nMatch++;
          break;
        }
      }
    }
  }

  if (do_cpu) {
    std::cout << nMatch << " seeds are matched" << std::endl;
    std::cout << "Matching rate: " << float(nMatch) / nSeed_cpu * 100 << "%"
              << std::endl;
  }

  if (!quiet) {
    if (do_cpu) {
      std::cout << "CPU Seed result:" << std::endl;

      for (auto& regionVec : seedVector_cpu) {
        for (size_t i = 0; i < regionVec.size(); i++) {
          const Acts::Seed<value_type>* seed = &regionVec[i];
          const value_type* sp = seed->sp()[0];
          std::cout << " (" << sp->x() << ", " << sp->y() << ", " << sp->z()
                    << ") ";
          sp = seed->sp()[1];
          std::cout << sp->sp()->surface << " (" << sp->x() << ", " << sp->y()
                    << ", " << sp->z() << ") ";
          sp = seed->sp()[2];
          std::cout << sp->sp()->surface << " (" << sp->x() << ", " << sp->y()
                    << ", " << sp->z() << ") ";
          std::cout << std::endl;
        }
      }

      std::cout << std::endl;
    }
    std::cout << "CUDA Seed result:" << std::endl;

    for (auto& regionVec : seedVector_cuda) {
      for (size_t i = 0; i < regionVec.size(); i++) {
        const Acts::Seed<value_type>* seed = &regionVec[i];
        const value_type* sp = seed->sp()[0];
        std::cout << " (" << sp->x() << ", " << sp->y() << ", " << sp->z()
                  << ") ";
        sp = seed->sp()[1];
        std::cout << sp->sp()->surface << " (" << sp->x() << ", " << sp->y()
                  << ", " << sp->z() << ") ";
        sp = seed->sp()[2];
        std::cout << sp->sp()->surface << " (" << sp->x() << ", " << sp->y()
                  << ", " << sp->z() << ") ";
        std::cout << std::endl;
      }
    }
  }

  std::cout << std::endl;
  std::cout << std::endl;

  return 0;
}
