// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MeasurementMapAlgorithm.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>

namespace ActsExamples {

MeasurementMapAlgorithm::MeasurementMapAlgorithm(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("MeasurementMapAlgorithm", std::move(logger)),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks");
  }
  if (m_cfg.outputMeasurementMap.empty()) {
    throw std::invalid_argument("Missing output measurement map");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputMeasurementMap.maybeInitialize(m_cfg.inputMeasurementMap);
  m_inputIndexRemapping.maybeInitialize(m_cfg.inputIndexRemapping);
  m_outputMeasurementMap.initialize(m_cfg.outputMeasurementMap);
}

ProcessCode MeasurementMapAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& tracks = m_inputTracks(ctx);

  // Start from the previous pass's map if provided, otherwise empty.
  UsedMeasurementMap outputMap;
  if (m_inputMeasurementMap.isInitialized()) {
    outputMap = m_inputMeasurementMap(ctx);
  }

  // Optional remapping from filtered-container indices to original indices.
  const MeasurementIndexRemapping* remapping = nullptr;
  if (m_inputIndexRemapping.isInitialized()) {
    remapping = &m_inputIndexRemapping(ctx);
  }

  std::size_t nUsed = 0;
  std::size_t nAlreadyPresent = 0;

  for (auto track : tracks) {
    for (auto state : track.trackStatesReversed()) {
      const bool isMeasurement = state.typeFlags().isMeasurement();
      const bool isOutlier =
          m_cfg.includeOutliers && state.typeFlags().isOutlier();
      if (!isMeasurement && !isOutlier) {
        continue;
      }
      if (!state.hasUncalibratedSourceLink()) {
        ACTS_WARNING("Measurement track state has no uncalibrated source link");
        continue;
      }

      Index filteredIdx =
          state.getUncalibratedSourceLink().get<IndexSourceLink>().index();

      // Translate to original-container index when a remapping is available.
      Index originalIdx =
          (remapping != nullptr) ? remapping->at(filteredIdx) : filteredIdx;

      auto [_, inserted] = outputMap.insert(originalIdx);
      if (inserted) {
        ++nUsed;
      } else {
        ++nAlreadyPresent;
      }
    }
  }

  ACTS_DEBUG("Recorded " << nUsed << " new measurement(s) from "
                         << tracks.size() << " track(s); " << nAlreadyPresent
                         << " already present from previous pass(es)");
  ACTS_DEBUG("Total used measurements: " << outputMap.size());

  m_outputMeasurementMap(ctx, std::move(outputMap));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
