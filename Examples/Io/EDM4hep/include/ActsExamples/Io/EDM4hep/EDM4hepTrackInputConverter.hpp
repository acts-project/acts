// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepInputConverter.hpp"

#include <string>

namespace podio {
class Frame;
}

namespace ActsExamples {

/// Read in a track collection as EDM4hep from a @c podio::Frame.
class EDM4hepTrackInputConverter : public EDM4hepInputConverter {
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

  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const final;

 private:
  Config m_cfg;

  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples
