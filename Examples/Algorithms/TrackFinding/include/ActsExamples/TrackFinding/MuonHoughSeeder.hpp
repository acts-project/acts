
#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Seeding/HoughTransformUtils.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/DriftCircle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "TH2D.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TCanvas.h"


namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples


namespace ActsExamples {
    
class MuonHoughSeeder final : public IAlgorithm {
 public:
  struct Config {
    std::string inSimHits;
    std::string inDriftCircles;
  };

  MuonHoughSeeder(Config cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;
  ProcessCode initialize() final;
  ProcessCode finalize() final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger& logger() const { return *m_logger; }


  ReadDataHandle<SimHitContainer> m_inputSimHits{this,
                                                           "InputSimHits"};
  ReadDataHandle<DriftCircleContainer> m_inputDriftCircles{this,
                                                           "InputDriftCircles"};
  std::unique_ptr<TCanvas> m_outCanvas;
  std::unique_ptr<TH2D> m_houghHist;
};

}  // namespace ActsExamples
