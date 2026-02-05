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

#include <memory>
#include <string>

namespace ActsPlugins::PodioUtil {
class ConversionHelper;
}

namespace podio {
class CollectionBase;
}

namespace ActsExamples {

/// Copy tracks from a ConstTrackContainer to a PodioTrackContainer
///
/// This algorithm reads a ConstTrackContainer from the whiteboard and
/// copies all tracks into a PodioTrackContainer, which is the Acts native
/// EDM4hep format. This preserves all track states and dynamic columns,
/// unlike the standard EDM4hep format conversion.
///
/// @note Currently a placeholder implementation pending resolution of
/// PodioTrackContainer collection release mechanism.
class TrackContainerToPodioConverter : public PodioOutputConverter {
 public:
  struct Config {
    /// Input track container
    std::string inputTracks;
    /// Output track collection in podio format
    std::string outputTracks = "tracks";
  };

  /// Constructor
  /// @param config is the configuration object
  /// @param helper is the conversion helper for surface/sourcelink mapping
  /// @param level is the output logging level
  explicit TrackContainerToPodioConverter(
      const Config& config,
      std::shared_ptr<ActsPlugins::PodioUtil::ConversionHelper> helper,
      Acts::Logging::Level level = Acts::Logging::INFO);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const final;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] context is the algorithm context for consistency
  ProcessCode execute(const AlgorithmContext& context) const final;

 private:
  Config m_cfg;

  std::shared_ptr<ActsPlugins::PodioUtil::ConversionHelper> m_helper;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};

  CollectionBaseWriteHandle m_outputTracks{this, "OutputTracks"};
  CollectionBaseWriteHandle m_outputTrackStates{this, "OutputTrackStates"};
};

}  // namespace ActsExamples
