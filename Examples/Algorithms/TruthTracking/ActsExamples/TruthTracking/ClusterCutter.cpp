#include "ActsExamples/TruthTracking/ClusterCutter.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace ActsExamples {

ClusterCutter::ClusterCutter(const Config& config, Acts::Logging::Level level)
    : IAlgorithm("ClusterCutter", level), m_cfg(config) {

  // Validate array bounds configuration
  if (m_cfg.maxR.size() != m_cfg.maxZ.size()) {
    throw std::invalid_argument("maxR and maxZ arrays must have the same length");
  }
  
  // Validate single value bounds
  if (m_cfg.minR < 0) {
    throw std::invalid_argument("minR must be non-negative");
  }
  
  // Validate input/output configuration
  if (m_cfg.inputClusters.empty()) {
    throw std::invalid_argument("Missing input clusters");
  }
  if (m_cfg.outputClusters.empty()) {
    throw std::invalid_argument("Missing output clusters");
  }

  m_inputClusters.initialize(m_cfg.inputClusters);
  m_outputClusters.initialize(m_cfg.outputClusters);
}

bool ClusterCutter::passesSelection(const Cluster& cluster) const {
  // Get cluster position (globalPosition is a member variable, not a function)
  const auto& globalPos = cluster.globalPosition;
  const auto r = std::sqrt(globalPos.x() * globalPos.x() + globalPos.y() * globalPos.y());
  const auto absZ = std::abs(globalPos.z());
  
  // Apply minimum bounds
  if (r < m_cfg.minR || absZ < m_cfg.minZ) {
    return false;
  }
  
  // Apply array-based bounds (similar to graph_modifier logic)
  bool passesArrayBounds = false;
  for (size_t i = 0; i < m_cfg.maxR.size(); ++i) {
    if (r <= m_cfg.maxR[i] && absZ <= m_cfg.maxZ[i]) {
      passesArrayBounds = true;
      break;
    }
  }
  
  return passesArrayBounds;
}

ProcessCode ClusterCutter::execute(const AlgorithmContext& ctx) const {
  const ClusterContainer& clusters = m_inputClusters(ctx);

  ACTS_VERBOSE("Total clusters available: " << clusters.size());

  ClusterContainer outputClusters;
  outputClusters.reserve(clusters.size());

  size_t selectedCount = 0;
  
  for (const auto& cluster : clusters) {
    if (passesSelection(cluster)) {
      outputClusters.push_back(cluster);
      selectedCount++;
      const auto& pos = cluster.globalPosition;
      ACTS_VERBOSE("Selected Cluster: " << pos.x() << " " << pos.y() << " " << pos.z() << " r=" << std::sqrt(pos.x()*pos.x() + pos.y()*pos.y()));
    } else {
      const auto& pos = cluster.globalPosition;
      ACTS_DEBUG("Removed Cluster: " << pos.x() << " " << pos.y() << " " << pos.z() << " r=" << std::sqrt(pos.x()*pos.x() + pos.y()*pos.y()));
    }
  }

  ACTS_DEBUG("Selected " << selectedCount << " out of " << clusters.size() << " clusters");

  m_outputClusters(ctx, std::move(outputClusters));

  return ProcessCode::SUCCESS;
}

} // namespace ActsExamples