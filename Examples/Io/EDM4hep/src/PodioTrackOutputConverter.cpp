// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/PodioTrackOutputConverter.hpp"

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "ActsPlugins/EDM4hep/PodioTrackContainer.hpp"
#include "ActsPlugins/EDM4hep/PodioTrackStateContainer.hpp"
#include "ActsPlugins/EDM4hep/PodioUtil.hpp"
#include "ActsPodioEdm/BoundParametersCollection.h"
#include "ActsPodioEdm/JacobianCollection.h"
#include "ActsPodioEdm/TrackCollection.h"
#include "ActsPodioEdm/TrackStateCollection.h"

#include <stdexcept>

namespace ActsExamples {

PodioTrackOutputConverter::PodioTrackOutputConverter(
    const Config& config,
    std::shared_ptr<ActsPlugins::PodioUtil::ConversionHelper> helper,
    Acts::Logging::Level level)
    : PodioOutputConverter("PodioTrackOutputConverter", level),
      m_cfg(config),
      m_helper(std::move(helper)) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output tracks collection");
  }
  if (!m_helper) {
    throw std::invalid_argument("Missing conversion helper");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
  m_outputTrackStates.initialize(m_cfg.outputTracks + "_trackStates");
}

ActsExamples::ProcessCode PodioTrackOutputConverter::execute(
    const AlgorithmContext& context) const {
  // Read input tracks
  const auto& inputTracks = m_inputTracks(context);

  ACTS_VERBOSE("Converting " << inputTracks.size()
                             << " tracks to Podio format");

  // Create external Podio collections that we'll own and write to the event
  // store
  auto trackCollection = std::make_unique<ActsPodioEdm::TrackCollection>();
  auto trackStateCollection =
      std::make_unique<ActsPodioEdm::TrackStateCollection>();
  auto paramsCollection =
      std::make_unique<ActsPodioEdm::BoundParametersCollection>();
  auto jacsCollection = std::make_unique<ActsPodioEdm::JacobianCollection>();

  // Create Podio backends using RefHolder for externally owned collections
  ActsPlugins::MutablePodioTrackStateContainer<Acts::RefHolder>
      trackStateContainer(*m_helper, *trackStateCollection, *paramsCollection,
                          *jacsCollection);
  ActsPlugins::MutablePodioTrackContainer<Acts::RefHolder> trackContainer(
      *m_helper, *trackCollection);

  // Wrap in Acts TrackContainer
  Acts::TrackContainer outputTracks{trackContainer, trackStateContainer};

  // Copy all tracks using the high-level copyFrom API
  for (const auto& inputTrack : inputTracks) {
    auto outputTrack = outputTracks.makeTrack();

    // Use copyFrom to handle all track-level and track state copying
    outputTrack.copyFrom(inputTrack);
  }

  ACTS_VERBOSE("Copied " << outputTracks.size() << " tracks with "
                         << trackStateContainer.size() << " track states");

  // Write collections to event store
  m_outputTracks(context, std::move(trackCollection));
  m_outputTrackStates(context, std::move(trackStateCollection));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> PodioTrackOutputConverter::collections() const {
  return {m_cfg.outputTracks, m_cfg.outputTracks + "_trackStates"};
}

}  // namespace ActsExamples
