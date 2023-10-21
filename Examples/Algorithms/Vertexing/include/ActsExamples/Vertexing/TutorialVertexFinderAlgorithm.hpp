// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <vector>

namespace Acts {
class MagneticFieldProvider;
}  // namespace Acts

namespace ActsExamples {
struct AlgorithmContext;

class TutorialVertexFinderAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Optional. Input track parameters collection
    std::string inputTrackParameters;
    /// Output proto vertex collection
    std::string outputProtoVertices;
    /// The magnetic field
    std::shared_ptr<Acts::MagneticFieldProvider> bField;
  };

  TutorialVertexFinderAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  /// Find vertices.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<TrackParametersContainer> m_inputTrackParameters{
      this, "InputTrackParameters"};

  WriteDataHandle<ProtoVertexContainer> m_outputProtoVertices{
      this, "OutputProtoVertices"};
};

}  // namespace ActsExamples
