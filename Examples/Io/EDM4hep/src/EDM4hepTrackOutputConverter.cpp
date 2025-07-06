// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepTrackOutputConverter.hpp"

#include "Acts/Plugins/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <stdexcept>

#include <edm4hep/TrackCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepTrackOutputConverter::EDM4hepTrackOutputConverter(
    const Config& config, Acts::Logging::Level level)
    : EDM4hepOutputConverter("EDM4hepTrackOutputConverter", level),
      m_cfg(config) {
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode EDM4hepTrackOutputConverter::execute(
    const AlgorithmContext& context) const {
  edm4hep::TrackCollection trackCollection;

  const auto& tracks = m_inputTracks(context);

  for (const auto& from : tracks) {
    auto to = trackCollection.create();
    Acts::EDM4hepUtil::writeTrack(context.geoContext, from, to, m_cfg.Bz);
  }

  m_outputTracks(context, std::move(trackCollection));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> EDM4hepTrackOutputConverter::collections() const {
  return {m_cfg.outputTracks};
}

}  // namespace ActsExamples
