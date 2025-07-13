// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Plugins/Root/RootMaterialTrackAccessor.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <TChain.h>

namespace ActsExamples {

/// @class RootMaterialTrackReader
///
/// @brief Reads in MaterialTrack information from a root file
/// and fills it into a format to be understood by the MaterialMapping
/// algorithm
class RootMaterialTrackReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// material collection to read
    std::string outputMaterialTracks = "material-tracks";
    /// name of the output tree
    std::string treeName = "material-tracks";
    /// List of input files
    std::vector<std::string> fileList;

    // Read surface information for the root file
    bool readCachedSurfaceInformation = false;
  };

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  RootMaterialTrackReader(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~RootMaterialTrackReader() override = default;

  /// Framework name() method
  std::string name() const override;

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const ActsExamples::AlgorithmContext& context) override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The logger
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_outputMaterialTracks{this, "OutputMaterialTracks"};

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// The number of events
  std::size_t m_events = 0;

  /// The batch size (number of track per events)
  std::size_t m_batchSize = 0;

  /// The input tree name
  std::unique_ptr<TChain> m_inputChain;

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> m_entryNumbers = {};

  Acts::RootMaterialTrackAccessor m_accessor;
};

}  // namespace ActsExamples
