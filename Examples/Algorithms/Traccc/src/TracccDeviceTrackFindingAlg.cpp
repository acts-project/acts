// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccDeviceTrackFindingAlg.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>

namespace ActsExamples {

TracccDeviceTrackFindingAlg::TracccDeviceTrackFindingAlg(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TracccDeviceTrackFindingAlg", std::move(logger)), m_cfg(cfg) {
  if (m_cfg.detector == nullptr) {
    throw std::invalid_argument(
        "TracccDeviceTrackFindingAlg: detector is not set");
  }
  if (m_cfg.magneticField == nullptr) {
    throw std::invalid_argument(
        "TracccDeviceTrackFindingAlg: magneticField is not set");
  }
  m_algorithm = m_cfg.detector->backend().makeCombinatorialKalmanFilter(
      m_cfg.finding, this->logger().clone("TracccTrackFinding"));

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode TracccDeviceTrackFindingAlg::execute(
    const AlgorithmContext& ctx) const {
  const TracccDetector& detector = *m_cfg.detector;

  TracccDeviceTracks tracks = (**m_algorithm)(
      detector.deviceDetector(), m_cfg.magneticField->deviceField(),
      m_inputMeasurements(ctx), m_inputTrackParameters(ctx));
  detector.backend().synchronize();

  ACTS_DEBUG("Found " << m_algorithm->copy().get_size(tracks.tracks)
                      << " tracks");

  m_outputTracks(ctx, std::move(tracks));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
