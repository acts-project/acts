// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/Seed.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/SpacePointData.hpp"
#include "Acts/Plugins/Sycl/Seeding/SeedFinder.hpp"
#include "Acts/Plugins/Sycl/Utilities/QueueWrapper.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/EventData/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>

#include <boost/type_erasure/any_cast.hpp>

#include "ATLASCuts.hpp"
#include "CommandLineArguments.h"
#include "SpacePoint.hpp"
#include "SpacePointContainer.hpp"
#include "vecmem/memory/sycl/device_memory_resource.hpp"
#include "vecmem/memory/sycl/host_memory_resource.hpp"

using namespace Acts::UnitLiterals;

auto readFile(const std::string& filename) -> std::vector<const SpacePoint*> {
  std::string line;
  std::vector<const SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()) {
    int id = 0;
    while (std::getline(spFile, line)) {
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      if (linetype == "lxyz") {
        float x;
        float y;
        float z;
        float r;
        float varianceR;
        float varianceZ;
        int layer;
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
        readSP.emplace_back(
            new SpacePoint(x, y, z, r, layer, varianceR, varianceZ, id));
        ++id;
      }
    }
  }
  return readSP;
}

template <typename external_spacepoint_t>
auto setupSeedFinderConfiguration()
    -> Acts::SeedFinderConfig<external_spacepoint_t> {
  Acts::SeedFinderConfig<external_spacepoint_t> config;
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
  // for sycl
  config.nTrplPerSpBLimit = 100;
  config.nAvgTrplPerSpBLimit = 6;
  return config;
}

auto setupSeedFinderOptions() {
  Acts::SeedFinderOptions options;
  options.bFieldInZ = 2_T;
  options.beamPos = {-.5_mm, -.5_mm};
  return options;
}

template <typename external_spacepoint_t>
auto setupSpacePointGridConfig(
    const Acts::SeedFinderConfig<external_spacepoint_t>& config,
    const Acts::SeedFinderOptions& options)
    -> std::pair<Acts::SpacePointGridConfig, Acts::SpacePointGridOptions> {
  Acts::SpacePointGridConfig gridConf{};
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;

  Acts::SpacePointGridOptions gridOpts{};
  gridOpts.bFieldInZ = options.bFieldInZ;
  return std::make_pair(gridConf, gridOpts);
}

auto main(int argc, char** argv) -> int {
  auto start_prep = std::chrono::system_clock::now();

  CommandLineArguments cmdlTool;
  cmdlTool.parse(argc, argv);

  if (!cmdlTool.inpFileName.empty() && !cmdlTool.inpFileExists) {
    std::cerr << "Input file not found\n";
    return -1;
  }

  auto spVec = readFile(cmdlTool.inpFileName);

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

  int numPhiNeighbors = 1;

  // extent used to store r range for middle spacepoint
  Acts::Extent rRangeSPExtent;

  const Acts::Range1D<float> rMiddleSPRange;

  std::vector<std::pair<int, int>> zBinNeighborsTop;
  std::vector<std::pair<int, int>> zBinNeighborsBottom;

  auto bottomBinFinder = std::make_unique<Acts::GridBinFinder<2ul>>(
      Acts::GridBinFinder<2ul>(numPhiNeighbors, zBinNeighborsBottom));
  auto topBinFinder = std::make_unique<Acts::GridBinFinder<2ul>>(
      Acts::GridBinFinder<2ul>(numPhiNeighbors, zBinNeighborsTop));
  auto config = setupSeedFinderConfiguration<value_type>();
  config = config.toInternalUnits().calculateDerivedQuantities();
  auto options = setupSeedFinderOptions();
  options = options.toInternalUnits().calculateDerivedQuantities(config);

  Acts::ATLASCuts<value_type> atlasCuts = Acts::ATLASCuts<value_type>();
  Acts::Sycl::DeviceExperimentCuts deviceAtlasCuts;
  config.seedFilter = std::make_unique<Acts::SeedFilter<value_type>>(
      Acts::SeedFilter<value_type>(Acts::SeedFilterConfig(), &atlasCuts));

  const Acts::Logging::Level logLvl =
      cmdlTool.csvFormat ? Acts::Logging::WARNING : Acts::Logging::INFO;
  Acts::Sycl::QueueWrapper queue(
      cmdlTool.deviceName,
      Acts::getDefaultLogger("Sycl::QueueWrapper", logLvl));
  vecmem::sycl::host_memory_resource resource(queue.getQueue());
  vecmem::sycl::device_memory_resource device_resource(queue.getQueue());
  Acts::Sycl::SeedFinder<value_type> syclSeedFinder(
      config, options, deviceAtlasCuts, queue, resource, &device_resource);

  Acts::SeedFinder<value_type> normalSeedFinder(config);
  auto [gridConfig, gridOpts] = setupSpacePointGridConfig(config, options);
  gridConfig = gridConfig.toInternalUnits();
  gridOpts = gridOpts.toInternalUnits();

  Acts::SpacePointGrid<value_type> grid =
      Acts::SpacePointGridCreator::createGrid<value_type>(gridConfig, gridOpts);
  Acts::SpacePointGridCreator::fillGrid(config, options, grid, spVec.begin(),
                                        spVec.end(), globalTool,
                                        rRangeSPExtent);

  std::array<std::vector<std::size_t>, 2ul> navigation;
  auto spGroup = Acts::BinnedSPGroup<value_type>(
      std::move(grid), *bottomBinFinder, *topBinFinder,
      std::move(navigation));

  auto end_prep = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsec_prep = end_prep - start_prep;
  double prepTime = elapsec_prep.count();

  if (!cmdlTool.csvFormat) {
    std::cout << "read " << spVec.size() << " SP from file "
              << cmdlTool.inpFileName << std::endl;
    std::cout << "Preparation time: " << std::to_string(prepTime) << std::endl;
  }

  // -------------------------------------- //
  // ----------- EXECUTE ON CPU ----------- //
  // -------------------------------------- //

  auto start_cpu = std::chrono::system_clock::now();
  uint group_count = 0;
  std::vector<std::vector<seed_type>> seedVector_cpu;

  if (!cmdlTool.onlyGpu) {
    decltype(normalSeedFinder)::SeedingState state;
    for (auto [bottom, middle, top] : spGroup) {
      normalSeedFinder.createSeedsForGroup(
          options, state, spGroup.grid(),
          std::back_inserter(seedVector_cpu.emplace_back()), bottom, middle,
          top, rMiddleSPRange);
      group_count++;
      if (!cmdlTool.allgroup && group_count >= cmdlTool.groups) {
        break;
      }
    }
  }

  auto end_cpu = std::chrono::system_clock::now();

  if (!cmdlTool.csvFormat) {
    std::cout << "Analyzed " << group_count << " groups for CPU" << std::endl;
  }

  // -------------------------------------- //
  // -------- EXECUTE ON GPU - SYCL ------- //
  // -------------------------------------- //

  auto start_sycl = std::chrono::system_clock::now();

  group_count = 0;
  std::vector<std::vector<seed_type>> seedVector_sycl;

  for (auto [bottom, middle, top] : spGroup) {
    seedVector_sycl.push_back(syclSeedFinder.createSeedsForGroup(
        spGroup.grid(), bottom, middle, top));
    group_count++;
    if (!cmdlTool.allgroup && group_count >= cmdlTool.groups) {
      break;
    }
  }
  auto end_sycl = std::chrono::system_clock::now();

  if (!cmdlTool.csvFormat) {
    std::cout << "Analyzed " << group_count << " groups for SYCL" << std::endl;
  }

  std::chrono::duration<double> elapsec_cpu = end_cpu - start_cpu;
  double cpuTime = elapsec_cpu.count();

  std::chrono::duration<double> elapsec_sycl = end_sycl - start_sycl;
  double syclTime = elapsec_sycl.count();

  auto textWidth = 20;
  auto numWidth = 11;

  int nSeed_cpu = 0;
  int nSeed_sycl = 0;
  int nMatch = 0;

  if (cmdlTool.matches && !cmdlTool.onlyGpu) {
    for (auto& outVec : seedVector_cpu) {
      nSeed_cpu += outVec.size();
    }

    for (auto& outVec : seedVector_sycl) {
      nSeed_sycl += outVec.size();
    }

    for (std::size_t i = 0; i < seedVector_cpu.size(); i++) {
      auto regionVec_cpu = seedVector_cpu[i];
      auto regionVec_sycl = seedVector_sycl[i];

      std::vector<std::vector<value_type>> seeds_cpu;
      std::vector<std::vector<value_type>> seeds_sycl;

      for (const auto& sd : regionVec_cpu) {
        std::vector<value_type> seed_cpu;
        seed_cpu.push_back(*(sd.sp()[0]));
        seed_cpu.push_back(*(sd.sp()[1]));
        seed_cpu.push_back(*(sd.sp()[2]));
        seeds_cpu.push_back(seed_cpu);
      }
      for (const auto& sd : regionVec_sycl) {
        std::vector<value_type> seed_sycl;
        seed_sycl.push_back(*(sd.sp()[0]));
        seed_sycl.push_back(*(sd.sp()[1]));
        seed_sycl.push_back(*(sd.sp()[2]));
        seeds_sycl.push_back(seed_sycl);
      }

      for (auto seed : seeds_cpu) {
        for (auto other : seeds_sycl) {
          if (*seed[0].sp() == *other[0].sp() and
              *seed[1].sp() == *other[1].sp() and
              *seed[2].sp() == *other[2].sp()) {
            nMatch++;
            break;
          }
        }
      }
    }
  }

  if (!cmdlTool.csvFormat) {
    std::cout << std::endl;
    std::cout
        << "------------------------- Time Metric -------------------------"
        << std::endl;
    std::cout << std::setw(textWidth) << " Device:";
    std::cout << std::setw(numWidth) << "CPU";
    std::cout << std::setw(numWidth) << "SYCL";
    std::cout << std::setw(textWidth) << "Speedup/ Agreement" << std::endl;
    std::cout << std::setw(textWidth) << " Time (s):";
    std::cout << std::setw(numWidth) << std::to_string(cpuTime);
    std::cout << std::setw(numWidth) << std::to_string(syclTime);
    std::cout << std::setw(textWidth) << std::to_string(cpuTime / syclTime);
    std::cout << std::endl;

    if (cmdlTool.matches && !cmdlTool.onlyGpu) {
      std::cout << std::setw(textWidth) << " Seeds found:";
      std::cout << std::setw(numWidth) << std::to_string(nSeed_cpu);
      std::cout << std::setw(numWidth) << std::to_string(nSeed_sycl);
      std::cout << std::setw(textWidth)
                << std::to_string(float(nMatch) / float(nSeed_cpu) * 100);
      std::cout << std::endl;
    }

    std::cout
        << "---------------------------------------------------------------"
        << std::endl;
    std::cout << std::endl;
  } else {
    std::cout << cpuTime << ',' << syclTime << ',' << cpuTime / syclTime << ','
              << nSeed_cpu << ',' << nSeed_sycl << ',' << nMatch << '\n';
  }

  for (const auto* S : spVec) {
    delete[] S;
  }

  return 0;
}
