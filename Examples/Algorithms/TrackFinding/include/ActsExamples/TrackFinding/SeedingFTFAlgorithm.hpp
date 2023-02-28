//basing off of SeedingOrtho and RosieTestFunc 

#pragma once

#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
//change these to FTF when done
#include "Acts/Seeding/SeedFinderFTF.hpp"
#include "Acts/Seeding/SeedFinderFTFConfig.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
   
#include <string>
#include <vector> 

namespace ActsExamples { 

class SeedingFTFAlgorithm final : public BareAlgorithm {
 public:
  struct Config {

    std::vector<std::string> inputSpacePoints;

    //two parameters that will be returned in function: 
    /// Output track seed collection.
    std::string outputSeeds;
    /// Output proto track collection.
    std::string outputProtoTracks;

    Acts::SeedFilterConfig seedFilterConfig;
    //change to FTF type
    Acts::SeedFinderFTFConfig<SimSpacePoint> seedFinderConfig;
    Acts::SeedFinderOptions seedFinderOptions;

  }; 


  //member functions of SFA class 


  //constructor: 
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingFTFAlgorithm(Config cfg, Acts::Logging::Level lvl); 
  
  //code to make the algorithm run 
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  //access to config 
  const Config& config() const { return m_cfg; }


 private: 
  //private memebers of SFA class 
  Config m_cfg; 
  Acts::SeedFinderFTF<SimSpacePoint> m_seedFinder; 
  //want to change to FTF version 

  //not copying print functions, think for debuging 
};

}  // namespace ActsExamples