// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MeasurementFilterAlgorithm.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>
#include <unordered_set>

namespace ActsExamples {

MeasurementFilterAlgorithm::MeasurementFilterAlgorithm(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("MeasurementFilterAlgorithm", std::move(logger)),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks");
  }
  if (m_cfg.inputMeasurementSubset.empty()) {
    throw std::invalid_argument("Missing input measurement subset");
  }
  if (m_cfg.outputMeasurementSubset.empty()) {
    throw std::invalid_argument("Missing output measurement subset");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputMeasurementSubset.initialize(m_cfg.inputMeasurementSubset);
  m_outputMeasurementSubset.initialize(m_cfg.outputMeasurementSubset);
}

ProcessCode MeasurementFilterAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& tracks = m_inputTracks(ctx);
  const auto& subset = m_inputMeasurementSubset(ctx);

  // Collect source-link indices consumed by this pass (always in
  // original-container space across all passes).
  std::unordered_set<Index> usedIndices;
  for (auto track : tracks) {
    for (auto state : track.trackStatesReversed()) {
      const bool isMeasurement = state.typeFlags().isMeasurement();
      const bool isOutlier =
          m_cfg.includeOutliers && state.typeFlags().isOutlier();
      if (!isMeasurement && !isOutlier) {
        continue;
      }
      if (!state.hasUncalibratedSourceLink()) {
        ACTS_WARNING("Track state has no uncalibrated source link");
        continue;
      }
      usedIndices.insert(
          state.getUncalibratedSourceLink().get<IndexSourceLink>().index());
    }
  }

  // Build output subset: collect indices not consumed in this pass.
  std::vector<MeasurementContainer::Index> validIndices;
  validIndices.reserve(subset.orderedIndices().size());

  std::size_t nFiltered = 0;
  for (const auto& sourceLink : subset.orderedIndices()) {
    if (usedIndices.count(sourceLink.index()) != 0u) {
      ++nFiltered;
      continue;
    }
    validIndices.push_back(sourceLink.index());
  }

  ACTS_DEBUG("Removed " << nFiltered << " measurement(s) used by "
                        << tracks.size() << " track(s); " << validIndices.size()
                        << " remaining");

  m_outputMeasurementSubset(
      ctx, MeasurementSubset(subset.container(), std::move(validIndices)));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
