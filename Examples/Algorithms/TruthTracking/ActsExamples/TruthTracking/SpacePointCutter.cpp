#include "ActsExamples/TruthTracking/SpacePointCutter.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace ActsExamples {

SpacePointCutter::SpacePointCutter(const Config& config, Acts::Logging::Level level)
    : IAlgorithm("SpacePointCutter", level), m_cfg(config) {

  // Validate array bounds configuration
  if (m_cfg.maxR.size() != m_cfg.maxZ.size()) {
    throw std::invalid_argument("maxR and maxZ arrays must have the same length");
  }
  
  // Validate single value bounds
  if (m_cfg.minR < 0) {
    throw std::invalid_argument("minR must be non-negative");
  }
  
  // Validate input/output configuration
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing input space points");
  }
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument("Missing output space points");
  }

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
}

bool SpacePointCutter::passesSelection(const SimSpacePoint& spacepoint) const {
  const auto r = spacepoint.r();
  const auto absZ = std::abs(spacepoint.z());
  
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

ProcessCode SpacePointCutter::execute(const AlgorithmContext& ctx) const {
  const SimSpacePointContainer& spacePoints = m_inputSpacePoints(ctx);

  ACTS_VERBOSE("Total space points available: " << spacePoints.size());

  SimSpacePointContainer outputSpacePoints;
  outputSpacePoints.reserve(spacePoints.size());

  size_t selectedCount = 0;
  
  for (const auto& spacePoint : spacePoints) {
    if (passesSelection(spacePoint)) {
      outputSpacePoints.push_back(spacePoint);
      selectedCount++;
      ACTS_VERBOSE("Selected SpacePoint: " << spacePoint.x() << " " << spacePoint.y() << " " << spacePoint.z() << " " << spacePoint.r());
    } else {
      ACTS_DEBUG("Removed SpacePoint: " << spacePoint.x() << " " << spacePoint.y() << " " << spacePoint.z() << " " << spacePoint.r());
    }
  }

  ACTS_DEBUG("Selected " << selectedCount << " out of " << spacePoints.size() << " space points");

  m_outputSpacePoints(ctx, std::move(outputSpacePoints));

  return ProcessCode::SUCCESS;
}

} // namespace ActsExamples