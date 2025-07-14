#pragma once

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/EventData/Cluster.hpp"

#include <string>
#include <limits>
#include <vector>

namespace ActsExamples {

class ClusterCutter final : public IAlgorithm {
public:
  struct Config {
    std::string inputClusters;
    std::string outputClusters;

    // Support both single values and arrays for R bounds
    std::vector<double> maxR = {std::numeric_limits<double>::infinity()};
    std::vector<double> maxZ = {std::numeric_limits<double>::infinity()};
    
    double minR = 0;
    double minZ = -std::numeric_limits<double>::infinity();
  };

  ClusterCutter(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

private:
  Config m_cfg;
  
  // Helper method to check if a cluster passes the selection criteria
  bool passesSelection(const Cluster& cluster) const;

  ReadDataHandle<ClusterContainer> m_inputClusters{this, "input_clusters"};
  WriteDataHandle<ClusterContainer> m_outputClusters{this, "output_clusters"};
};

} // namespace ActsExamples