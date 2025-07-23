#include "ActsExamples/Digitization/MuonSpacePointDigitizer.hpp"
    
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

namespace ActsExamples {
MuonSpacePointDigitizer::MuonSpacePointDigitizer(
    const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("MuonSpacePointDigitizer", lvl),
      m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("MuonSpacePointDigitizer", lvl)) {}


ProcessCode MuonSpacePointDigitizer::initialize() {
    return ProcessCode::SUCCESS;
}

ProcessCode MuonSpacePointDigitizer::execute(
    const AlgorithmContext& ctx) const {
  return ProcessCode::SUCCESS;
}

}