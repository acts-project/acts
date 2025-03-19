// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepTrackInputConverter.hpp"

#include "Acts/Plugins/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <stdexcept>

#include <edm4hep/TrackCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepTrackInputConverter::EDM4hepTrackInputConverter(
    const Config& config, Acts::Logging::Level level)
    : EDM4hepInputConverter("EDM4hepTrackInputConverter", level,
                            config.inputFrame),
      m_cfg(config) {
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode EDM4hepTrackInputConverter::convert(
    const AlgorithmContext& ctx, const podio::Frame& frame) const {
  const auto& trackCollection =
      frame.get<edm4hep::TrackCollection>(m_cfg.inputTracks);

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  for (const auto& inputTrack : trackCollection) {
    auto track = tracks.makeTrack();
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
