// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioCollectionDataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"
#include "ActsPodioEdm/TrackCollection.h"
#include "ActsPodioEdm/TrackStateCollection.h"
#include "ActsPodioEdm/TrackStateHitLinkCollection.h"
#include "ActsPodioEdm/TrackerHitLocalCollection.h"

#include <memory>
#include <string>

namespace ActsExamples {

/// Write out a track container to PODIO track collections
///
/// This algorithm reads a ConstTrackContainer from the whiteboard and
/// writes it to PodioTrackContainer format, which is the Acts native
/// EDM4hep format. This preserves all track states and dynamic columns,
/// unlike the standard EDM4hep format conversion.
///
/// The TrackerHitLocal collection must be written first by
/// EDM4hepMeasurementOutputConverter and supplied via inputTrackerHitsLocal.
class PodioTrackOutputConverter : public PodioOutputConverter {
 public:
  struct Config {
    /// Input track container
    std::string inputTracks;
    /// Output track collection in podio format
    std::string outputTracks;
    /// Input TrackerHitLocal collection
    std::string inputTrackerHitsLocal;
    /// DD4hep detector
    std::shared_ptr<DD4hepDetector> detector;
  };

  /// Constructor
  /// @param config is the configuration object
  /// @param logger is the logger
  explicit PodioTrackOutputConverter(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const final;

 protected:
  ProcessCode execute(const AlgorithmContext& context) const final;

 private:
  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  PodioCollectionReadHandle<ActsPodioEdm::TrackerHitLocalCollection>
      m_inputTrackerHitsLocal{this, "InputTrackerHitsLocal"};

  PodioCollectionWriteHandle<ActsPodioEdm::TrackCollection> m_outputTracks{
      this, "OutputTracks"};
  PodioCollectionWriteHandle<ActsPodioEdm::TrackStateCollection>
      m_outputTrackStates{this, "OutputTrackStates"};
  PodioCollectionWriteHandle<ActsPodioEdm::TrackStateHitLinkCollection>
      m_outputTrackStateHitLinks{this, "OutputTrackStateHitLinks"};
};

}  // namespace ActsExamples
