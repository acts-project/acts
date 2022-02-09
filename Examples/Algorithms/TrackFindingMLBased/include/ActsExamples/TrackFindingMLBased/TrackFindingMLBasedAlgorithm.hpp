#pragma once

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFinding.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class TrackFindingMLBasedAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input spacepoints collection.
    std::string inputSpacePoints;

    /// Output protoTracks collection.
    std::string outputProtoTracks;

    /// ML based track finder
    std::shared_ptr<Acts::ExaTrkXTrackFinding> trackFinderML;

    // NOTE the other config parameters for the Exa.TrkX class for now are just
    // initialized as the defaults
  };

  /// Constructor of the track finding algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  TrackFindingMLBasedAlgorithm(Config cfg, Acts::Logging::Level lvl);

  virtual ~TrackFindingMLBasedAlgorithm() {}

  /// Framework execute method of the track finding algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  // configuration
  Config m_cfg;
};

}  // namespace ActsExamples
