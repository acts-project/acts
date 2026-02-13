// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/EventData/Vertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace ActsExamples {

/// Matches reconstructed vertices to truth interactions and vice versa
class VertexTruthMatcher final : public IAlgorithm {
 public:
  struct Config {
    /// Input vertex collection.
    std::string inputVertices;
    /// Tracks object from track finding.
    std::string inputTracks;
    /// All input truth particle collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Output vertex-truth matching.
    std::string outputVertexTruthMatching;
    /// Output truth-vertex matching.
    std::string outputTruthVertexMatching;

    /// Minimum fraction of track weight matched between truth
    /// and reco vertices to consider as truth matched.
    double vertexMatchThreshold = 0.7;
    /// Minimum track weight for track to be considered as part of the fit
    double minTrkWeight = 0.1;
  };

  VertexTruthMatcher(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<VertexContainer> m_inputVertices{this, "InputVertices"};
  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};

  WriteDataHandle<std::vector<VertexToTruthMatching>>
      m_outputVertexTruthMatching{this, "OutputVertexTruthMatching"};
  WriteDataHandle<std::map<SimVertexBarcode, VertexToRecoMatching>>
      m_outputTruthVertexMatching{this, "OutputTruthVertexMatching"};
};

}  // namespace ActsExamples
