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
  try {
    po::options_description optionsDescription("Allowed options");
    optionsDescription.add_options()
      ("help,h", "Print usage message.")
      ("inputfile,f", po::value<std::string>(), "Provide path for input file.")
      ("platforms,p","List available platforms and devices.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, optionsDescription), vm);

    if (vm.count("help") != 0) {  
      std::cout << optionsDescription << "\n";
      return 0;
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

      for (; !(groupIt == spGroup.end()); ++groupIt) {
        seedVector_cpu.push_back(normalSeedfinder.createSeedsForGroup(
            groupIt.bottom(), groupIt.middle(), groupIt.top()));
        group_count++;
        if (!allgroup && group_count >= 500) {
          break;
        }
      }

      std::cout << "Analyzed " << group_count << " groups for CPU" << std::endl;

      auto end_cpu = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsec_cpu = end_cpu - start_cpu;
      double cpuTime = elapsec_cpu.count();

      //----------- SYCL ----------//

      auto start_sycl = std::chrono::system_clock::now();

      group_count = 0;
      std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector_sycl;
      groupIt = spGroup.begin();

      for (; !(groupIt == spGroup.end()); ++groupIt) {
        seedVector_sycl.push_back(syclSeedfinder.createSeedsForGroup(
            groupIt.bottom(), groupIt.middle(), groupIt.top()));
        group_count++;
        if (!allgroup && group_count >= 500){
            break;
        }
      }
      auto end_sycl = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsec_sycl = end_sycl - start_sycl;
      double syclTime = elapsec_sycl.count();

      std::cout << "Analyzed " << group_count << " groups for CUDA" << std::endl;

      std::cout << std::endl;
      std::cout << "----------------------- Time Metric -----------------------" << std::endl;
      std::cout << "Seedfinding_Time  " << std::setw(11)
                << (std::to_string(cpuTime)) << "  " << std::setw(11);
      std::cout << std::to_string(syclTime) << std::endl;
      std::cout << "-----------------------------------------------------------"
                << std::endl;
      std::cout << std::endl;

    for(const auto *S: spVec) {
        delete[] S;
      }
    }

    if(vm.count("platforms") != 0){
      Acts::Sycl::testDevice();
      // Acts::Sycl::outputPlatforms();
    }
  }
  catch (std::exception &e){
    std::cerr << e.what() << std::endl;
  }

  return 0;
}


