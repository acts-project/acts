#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFinding.hpp"

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

      ExaTrkXTrackFinding::Config exaTrkxConfig;

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

    ExaTrkXTrackFinding m_exaTrkx;

   private:
    // configuration
    Config m_cfg;
};

}