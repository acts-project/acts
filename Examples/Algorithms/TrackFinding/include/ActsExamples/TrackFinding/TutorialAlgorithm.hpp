#pragma once

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>

namespace ActsExamples {

class TutorialAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// The message we are going to log.
    std::string message = "hello world!";
  };

  TutorialAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  /// This function will be called on each event by the sequencer.
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
