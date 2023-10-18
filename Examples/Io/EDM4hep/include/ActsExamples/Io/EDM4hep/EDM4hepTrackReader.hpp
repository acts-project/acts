// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/SequenceElement.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <string>

#include <podio/ROOTFrameReader.h>

namespace ActsExamples {

class EDM4hepTrackReader : public IReader {
 public:
  struct Config {
    /// Input file path
    std::string inputPath;
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
  EDM4hepTrackReader(const Config& config,
                     Acts::Logging::Level level = Acts::Logging::INFO);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  std::string name() const final;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final;

  /// Read out data from the input stream.
  ProcessCode read(const ActsExamples::AlgorithmContext& ctx) final;

 private:
  Config m_cfg;

  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};

  podio::ROOTFrameReader m_reader;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
