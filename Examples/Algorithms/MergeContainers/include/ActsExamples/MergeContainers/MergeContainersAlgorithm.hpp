// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ActsExamples {
struct AlgorithmContext;

class MergeContainersAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input track parameters collection
    std::vector<std::string>  inputTrackParameters;
    /// Input track parameters collection
    std::string outputTrackParameters;
    /// Input track parameters collection
    std::vector<std::string>  inputTracks;
    /// Input track parameters collection
    std::string outputTracks;
  };

  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  MergeContainersAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the ambiguity resolution algorithm.
  ///
  /// @param cxt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }
  

 private:
  Config m_cfg;

  std::vector<std::unique_ptr<ReadDataHandle<TrackParametersContainer>>> m_inputTrackParameters{};
  WriteDataHandle<TrackParametersContainer> m_outputTrackParameters{this, "OutputTrackParameters"};

  std::vector<std::unique_ptr<ReadDataHandle<ConstTrackContainer>>> m_inputTracks{};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};

};

}  // namespace ActsExamples
