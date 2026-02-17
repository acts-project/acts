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
#include "ActsExamples/Io/Podio/CollectionBaseWriteHandle.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"

#include <string>

namespace podio {
class CollectionBase;
}

namespace ActsExamples {

/// Write out a track collection to EDM4hep objects
class EDM4hepTrackOutputConverter : public PodioOutputConverter {
 public:
  struct Config {
    /// Input track collection
    std::string inputTracks;
    /// Output track collection in edm4hep
    std::string outputTracks = "ActsTracks";
    /// Magnetic field along the z axis (needed for the conversion of
    /// parameters)
    double Bz{};
  };

  /// constructor
  /// @param config is the configuration object
  /// @param level is the output logging level
  explicit EDM4hepTrackOutputConverter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const final;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] context is the algorithm context for consistency
  ProcessCode execute(const AlgorithmContext& context) const final;

 private:
  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};

  CollectionBaseWriteHandle m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples
