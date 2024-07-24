// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/FilterMeasurements/FilterMeasurementsAlgorithm.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <unordered_map>

ActsExamples::FilterMeasurementsAlgorithm::
    FilterMeasurementsAlgorithm(
        ActsExamples::FilterMeasurementsAlgorithm::Config cfg,
        Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("FilterMeasurementsAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
        
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links input collection");
  }  
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing tracks input collection");
  }
  if (m_cfg.outputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links output collection");
  }

  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
}

ActsExamples::ProcessCode
ActsExamples::FilterMeasurementsAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& tracks = m_inputTracks(ctx);
  IndexSourceLinkContainer outputSourceLinks = m_inputSourceLinks(ctx);

  ACTS_VERBOSE("Number of input tracks: " << tracks.size());


  for (const auto& track : tracks) {
    for (auto ts : track.trackStatesReversed()) {
      if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        Acts::SourceLink sourceLink = ts.getUncalibratedSourceLink();
        const auto islink = sourceLink.get<IndexSourceLink>();
        outputSourceLinks.erase(std::remove(outputSourceLinks.begin(),outputSourceLinks.end(), islink), outputSourceLinks.end());
      }
    }
  }

  m_outputSourceLinks(ctx, std::move(outputSourceLinks));
  return ActsExamples::ProcessCode::SUCCESS;

}