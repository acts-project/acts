// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/SequenceElement.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <string>

namespace podio {
class Frame;
}

namespace ActsExamples {

class EDM4hepTrackInputConverter : public IAlgorithm {
 public:
  struct Config {
    std::string inputFrame;
    /// Input track collection name in edm4hep
    std::string inputTracks = "ActsTracks";
    /// Output track collection
    std::string outputTracks;
    /// Magnetic field along the z axis (needed for the conversion of
    /// parameters)
    double Bz;
  };

  /// constructor
  /// @param config is the configuration object
  /// @param level is the output logging level
  explicit EDM4hepTrackInputConverter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
  Config m_cfg;

  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
  ReadDataHandle<podio::Frame> m_inputFrame{this, "InputFrame"};
};

}  // namespace ActsExamples
