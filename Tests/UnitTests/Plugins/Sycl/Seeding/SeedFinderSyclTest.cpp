#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <memory>
#include <cmath>

#include <boost/type_erasure/any_cast.hpp>
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include "ATLASCuts.hpp"
#include "SpacePoint.hpp"
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"

// use boost program options for parsing command line arguments
#include <boost/program_options.hpp>
#include <string>
namespace po = boost::program_options;

auto readFile(const std::string& filename) -> std::vector<const SpacePoint*> {
  std::string line;
  std::vector<const SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()) {
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
        readSP.emplace_back(new SpacePoint(x, y, z, r, layer, varianceR, varianceZ));
      }
    }
  }
  return readSP;
}

template <typename external_spacepoint_t>
auto setupSeedfinderConfiguration() -> Acts::SeedfinderConfig<external_spacepoint_t> {
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

  // for sycl
  config.nTrplPerSpBLimit = 100;
  config.nAvgTrplPerSpBLimit = 6;
  return config;
}

template <typename external_spacepoint_t>
auto setupSpacePointGridConfig(const Acts::SeedfinderConfig<external_spacepoint_t> &config) -> Acts::SpacePointGridConfig {
  Acts::SpacePointGridConfig gridConf{};
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;
  return gridConf;
}

auto main(int argc, char** argv) -> int {
  bool allgroup(false);
  bool cpu(true);
  bool matches(false);
  int groups(500);
  try {
    po::options_description optionsDescription("Allowed options");
    optionsDescription.add_options()
      ("help,h", "Print usage message.")
      ("inputfile,f", po::value<std::string>(), "Provide path for input file.")
      ("platforms,p","List available platforms and devices.")
      ("no_cpu,c","Do not execute code on cpu")
      ("groups,g",po::value<int>(),"Add number of groups to execute on")
      ("matches,m","Count matches")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, optionsDescription), vm);

    if (vm.count("help") != 0) {  
      std::cout << optionsDescription << "\n";
      return 0;
    }

    if(vm.count("no_cpu") != 0) {  
      cpu = false;
    }

    if(vm.count("matches") != 0) {
      matches = true;
    }

    if(vm.count("groups") != 0){
      groups = vm["groups"].as<int>();
    }

    if(vm.count("inputfile") != 0){
      std::string filename = vm["inputfile"].as<std::string>();
      auto spVec = readFile(filename);

      auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
        Acts::BinFinder<SpacePoint>());
      auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
        Acts::BinFinder<SpacePoint>());
      auto config = setupSeedfinderConfiguration<SpacePoint>();
      Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
      config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
        Acts::SeedFilter<SpacePoint>(Acts::SeedFilterConfig(), &atlasCuts));
      Acts::Sycl::Seedfinder<SpacePoint> syclSeedfinder(config);
      Acts::Seedfinder<SpacePoint> normalSeedfinder(config);
      auto covarianceTool = [=](const SpacePoint& sp, float /*unused*/, float /*unused*/, float /*unused*/) -> Acts::Vector2D {
        return {sp.varianceR, sp.varianceZ};
      };
      std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid = 
        Acts::SpacePointGridCreator::createGrid<SpacePoint>(setupSpacePointGridConfig(config));

      auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(), spVec.end(), covarianceTool,
                                                bottomBinFinder, topBinFinder, std::move(grid), config);
      std::cout << "read " << spVec.size() << " SP from file " << filename << std::endl;

      // ********* EXECUTE ON CPU ********** //

      auto start_cpu = std::chrono::system_clock::now();
      int group_count = 0;
      auto groupIt = spGroup.begin();

      std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector_cpu;
      if(cpu) {

        for (; !(groupIt == spGroup.end()); ++groupIt) {
          seedVector_cpu.push_back(normalSeedfinder.createSeedsForGroup(
              groupIt.bottom(), groupIt.middle(), groupIt.top()));
          group_count++;
          if (!allgroup && group_count >= groups) {
            break;
          }
        }
      }

      std::cout << "Analyzed " << group_count << " groups for CPU" << std::endl;

      auto end_cpu = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsec_cpu = end_cpu - start_cpu;
      double cpuTime = elapsec_cpu.count();

      //----------- EXECUTE ON GPU - SYCL ----------//

      auto start_sycl = std::chrono::system_clock::now();

      group_count = 0;
      std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector_sycl;
      groupIt = spGroup.begin();

      for (; !(groupIt == spGroup.end()); ++groupIt) {
        seedVector_sycl.push_back(syclSeedfinder.createSeedsForGroup(
            groupIt.bottom(), groupIt.middle(), groupIt.top()));
        group_count++;
        if (!allgroup && group_count >= groups){
            break;
        }
      }
      auto end_sycl = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsec_sycl = end_sycl - start_sycl;
      double syclTime = elapsec_sycl.count();

      std::cout << "Analyzed " << group_count << " groups for SYCL" << std::endl;

      std::cout << std::endl;
      std::cout << "----------------------- Time Metric -----------------------" << std::endl;
      std::cout << std::setw(20) << " Device:" << std::setw(11) << "CPU";
      std::cout << std::setw(11) << "SYCL";
      std::cout << std::setw(11) << "speedup" << std::endl;
      std::cout << std::setw(20) << " Seedfinding_Time:";
      std::cout << std::setw(11) << std::to_string(cpuTime) << " ";
      std::cout << std::setw(11) << std::to_string(syclTime);
      std::cout << std::setw(11) << std::to_string(cpuTime/syclTime);
      std::cout << std::endl;
      std::cout << "-----------------------------------------------------------" << std::endl;

      for(const auto *S: spVec) {
        delete[] S;
      }

      if(matches) {
        int nSeed_cpu = 0;
        for (auto& outVec : seedVector_cpu) {
          nSeed_cpu += outVec.size();
        }

        int nSeed_cuda = 0;
        for (auto& outVec : seedVector_sycl) {
          nSeed_cuda += outVec.size();
        }

        std::cout << "Number of Seeds (CPU | CUDA): " << nSeed_cpu << " | "
                  << nSeed_cuda << std::endl;

        int nMatch = 0;

        for (size_t i = 0; i < seedVector_cpu.size(); i++) {
          auto regionVec_cpu = seedVector_cpu[i];
          auto regionVec_cuda = seedVector_sycl[i];

          std::vector<std::vector<SpacePoint>> seeds_cpu;
          std::vector<std::vector<SpacePoint>> seeds_cuda;

          // for (size_t i_cpu = 0; i_cpu < regionVec_cpu.size(); i_cpu++) {
          for (auto sd : regionVec_cpu) {
            std::vector<SpacePoint> seed_cpu;
            seed_cpu.push_back(*(sd.sp()[0]));
            seed_cpu.push_back(*(sd.sp()[1]));
            seed_cpu.push_back(*(sd.sp()[2]));

            seeds_cpu.push_back(seed_cpu);
          }

          for (auto sd : regionVec_cuda) {
            std::vector<SpacePoint> seed_cuda;
            seed_cuda.push_back(*(sd.sp()[0]));
            seed_cuda.push_back(*(sd.sp()[1]));
            seed_cuda.push_back(*(sd.sp()[2]));

            seeds_cuda.push_back(seed_cuda);
          }

          for (auto seed : seeds_cpu) {
            for (auto other : seeds_cuda) {
              if (seed[0] == other[0] && seed[1] == other[1] && seed[2] == other[2]) {
                nMatch++;
                break;
              }
            }
          }
        }

        if (cpu) {
          std::cout << nMatch << " seeds are matched" << std::endl;
          std::cout << "Matching rate: " << float(nMatch) / nSeed_cpu * 100 << "%"
                    << std::endl;
        }
      }
    }

    if(vm.count("platforms") != 0){
      Acts::Sycl::testDevice();
      Acts::Sycl::outputPlatforms();
    }
  }
  catch (std::exception &e){
    std::cerr << e.what() << std::endl;
  }

  return 0;
}


