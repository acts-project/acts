// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepTrackReader.hpp"

#include "Acts/Plugins/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <stdexcept>

#include <edm4hep/TrackCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepTrackReader::EDM4hepTrackReader(const Config& config,
                                       Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepSimHitReader", level)) {
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output trajectories collection");
  }

  m_outputTracks.initialize(m_cfg.outputTracks);

  m_reader.openFile(m_cfg.inputPath);
}

std::pair<size_t, size_t> EDM4hepTrackReader::availableEvents() const {
  return {0, m_reader.getEntries("events")};
}

std::string EDM4hepTrackReader::EDM4hepTrackReader::name() const {
  return "EDM4hepTrackReader";
}

ProcessCode EDM4hepTrackReader::read(const AlgorithmContext& ctx) {
  podio::Frame frame = m_reader.readEntry("events", ctx.eventNumber);

  const auto& trackCollection =
      frame.get<edm4hep::TrackCollection>(m_cfg.inputTracks);

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  for (const auto& inputTrack : trackCollection) {
    auto track = tracks.getTrack(tracks.addTrack());
    Acts::EDM4hepUtil::readTrack(inputTrack, track, m_cfg.Bz);
  }

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};

  m_outputTracks(ctx, std::move(constTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
