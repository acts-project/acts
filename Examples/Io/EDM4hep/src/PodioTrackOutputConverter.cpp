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
#include "ActsExamples/Io/EDM4hep/DD4hepPodioConversionHelper.hpp"
#include "ActsPlugins/EDM4hep/PodioTrackContainer.hpp"
#include "ActsPlugins/EDM4hep/PodioTrackStateContainer.hpp"
#include "ActsPlugins/EDM4hep/PodioUtil.hpp"
#include "ActsPodioEdm/BoundParametersCollection.h"
#include "ActsPodioEdm/JacobianCollection.h"
#include "ActsPodioEdm/TrackCollection.h"
#include "ActsPodioEdm/TrackStateCollection.h"
#include "ActsPodioEdm/TrackStateHitLinkCollection.h"

#include <stdexcept>

namespace ActsExamples {

PodioTrackOutputConverter::PodioTrackOutputConverter(
    const Config& config, std::unique_ptr<const Acts::Logger> logger)
    : PodioOutputConverter("PodioTrackOutputConverter", std::move(logger)),
      m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output tracks collection");
  }
  if (m_cfg.inputTrackerHitsLocal.empty()) {
    throw std::invalid_argument("Missing input tracker hits local collection");
  }
  if (m_cfg.detector == nullptr) {
    throw std::invalid_argument("Missing detector");
  }
  if (m_cfg.detector->trackingGeometry() == nullptr) {
    throw std::invalid_argument("Tracking geometry not found");
  }

  const auto prefix = m_cfg.outputTracks;

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputTrackerHitsLocal.initialize(m_cfg.inputTrackerHitsLocal);
  m_outputTracks.initialize(prefix);
  m_outputTrackStates.initialize(
      ActsPlugins::PodioTrackStateContainerBase::trackStatesKey(prefix));
  m_outputTrackStateHitLinks.initialize(
      ActsPlugins::PodioTrackStateContainerBase::trackStateHitLinksKey(prefix));
}

ActsExamples::ProcessCode PodioTrackOutputConverter::execute(
    const AlgorithmContext& context) const {
  const auto& inputTracks = m_inputTracks(context);

  ACTS_VERBOSE("Converting " << inputTracks.size()
                             << " tracks to Podio format");

  const auto& trackerHitsLocal = m_inputTrackerHitsLocal(context);

  auto trackStateHitLinks =
      std::make_unique<ActsPodioEdm::TrackStateHitLinkCollection>();

  DD4hepPodioConversionHelper helper(*m_cfg.detector, trackerHitsLocal);

  auto trackCollection = std::make_unique<ActsPodioEdm::TrackCollection>();
  auto trackStateCollection =
      std::make_unique<ActsPodioEdm::TrackStateCollection>();
  auto paramsCollection =
      std::make_unique<ActsPodioEdm::BoundParametersCollection>();
  auto jacsCollection = std::make_unique<ActsPodioEdm::JacobianCollection>();

  ActsPlugins::MutablePodioTrackStateContainer trackStateContainer(
      helper, *trackStateCollection, *paramsCollection, *jacsCollection,
      trackStateHitLinks.get());
  ActsPlugins::MutablePodioTrackContainer trackContainer(helper,
                                                         *trackCollection);

  Acts::TrackContainer outputTracks{trackContainer, trackStateContainer};

  std::size_t nUncalibratedSourceLinks = 0;

  for (const auto& inputTrack : inputTracks) {
    for (const auto& ts : inputTrack.trackStatesReversed()) {
      if (ts.hasUncalibratedSourceLink()) {
        nUncalibratedSourceLinks++;
      }
    }

    auto outputTrack = outputTracks.makeTrack();
    outputTrack.copyFrom(inputTrack);
  }

  ACTS_VERBOSE("Copied " << outputTracks.size() << " tracks with "
                         << trackStateContainer.size() << " track states");

  std::size_t nUncalibratedSourceLinksInOutput = 0;
  for (const auto& outputTrack : outputTracks) {
    for (const auto& ts : outputTrack.trackStates()) {
      if (ts.hasUncalibratedSourceLink()) {
        nUncalibratedSourceLinksInOutput++;
      }
    }
  }

  if (nUncalibratedSourceLinks != nUncalibratedSourceLinksInOutput) {
    ACTS_ERROR(
        "Number of uncalibrated source links in input and output do not match: "
        << nUncalibratedSourceLinks
        << " != " << nUncalibratedSourceLinksInOutput);
    return ProcessCode::ABORT;
  }

  m_outputTracks(context, std::move(trackCollection));
  m_outputTrackStates(context, std::move(trackStateCollection));
  m_outputTrackStateHitLinks(context, std::move(trackStateHitLinks));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> PodioTrackOutputConverter::collections() const {
  const auto prefix = m_cfg.outputTracks;
  return {
      prefix,
      ActsPlugins::PodioTrackStateContainerBase::trackStatesKey(prefix),
      ActsPlugins::PodioTrackStateContainerBase::trackStateHitLinksKey(prefix),
  };
}

}  // namespace ActsExamples
