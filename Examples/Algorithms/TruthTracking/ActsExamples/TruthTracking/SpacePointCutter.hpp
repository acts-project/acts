#pragma once

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <string>
#include <limits>
#include <vector>

namespace ActsExamples {

class SpacePointCutter final : public IAlgorithm {
public:
  struct Config {
    std::string inputSpacePoints;
    std::string outputSpacePoints;

    // Support both single values and arrays for R bounds
    std::vector<double> maxR = {std::numeric_limits<double>::infinity()};
    std::vector<double> maxZ = {std::numeric_limits<double>::infinity()};
    
    double minR = 0;
    double minZ = -std::numeric_limits<double>::infinity();
  };

  SpacePointCutter(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

private:
  Config m_cfg;
  
  // Helper method to check if a point passes the selection criteria
  bool passesSelection(const SimSpacePoint& spacepoint) const;

  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this, "input_space_points"};
  WriteDataHandle<SimSpacePointContainer> m_outputSpacePoints{this, "output_space_points"};
};

} // namespace ActsExamples
