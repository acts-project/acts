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
#include "Acts/Seeding/SpacePointGrid.hpp"

#include "ATLASCuts.hpp"
#include "SpacePoint.hpp"
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"

// use boost program options for parsing command line arguments
#include <boost/program_options.hpp>
namespace po = boost::program_options;

auto readFile(const std::string& filename) -> std::vector<std::unique_ptr<SpacePoint>> {
  std::string line;
  std::vector<std::unique_ptr<SpacePoint>> readSP;

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
        readSP.emplace_back(std::make_unique<SpacePoint>(x, y, z, r, layer, varianceR, varianceZ));
      }
    }
  }
  return readSP;
}

template <typename external_spacepoint_t>
auto instatiateSeedfinderConfiguration() -> Acts::SeedfinderConfig<external_spacepoint_t> {
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

auto main(int argc, char** argv) -> int {
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

      auto config = instatiateSeedfinderConfiguration<SpacePoint>();
      Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
      config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
        Acts::SeedFilter<SpacePoint>(Acts::SeedFilterConfig(), &atlasCuts));
      auto SF(config);

      auto spVec = readFile(filename);
      std::cout << "read " << spVec.size() << " SP from file " << filename << std::endl;
    }

    if(vm.count("platforms") != 0){
      Acts::Sycl::outputPlatforms();
    }
  }
  catch (std::exception &e){
    std::cerr << e.what() << std::endl;
  }

  return 0;
}