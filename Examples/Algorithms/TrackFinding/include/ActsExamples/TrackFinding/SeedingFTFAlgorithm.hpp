//basing off of SeedingOrtho and RosieTestFunc 

#pragma once

#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
//change these to FTF when done
#include "Acts/Seeding/SeedFinderFTF.hpp"
#include "Acts/Seeding/SeedFinderFTFConfig.hpp"
//in core 
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
   
#include <string>
#include <vector> 

namespace ActsExamples { 

class SeedingFTFAlgorithm final : public IAlgorithm {
 public:
  struct Config {

    std::vector<std::string> inputSpacePoints;

    //two parameters that will be returned in function: 
    /// Output track seed collection.
    std::string outputSeeds;

    Acts::SeedFilterConfig seedFilterConfig;
    //change to FTF type
    Acts::SeedFinderFTFConfig<SimSpacePoint> seedFinderConfig;
    Acts::SeedFinderOptions seedFinderOptions;

    std::string layerMappingFile ; 

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

  //new since update 

  std::vector<std::unique_ptr<ReadDataHandle<SimSpacePointContainer>>>
      m_inputSpacePoints{};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};

  //not copying print functions, think for debuging 
};

}  // namespace ActsExamples