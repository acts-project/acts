// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccDeviceTrackParamEstimationAlg.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>

namespace ActsExamples {

TracccDeviceTrackParamEstimationAlg::TracccDeviceTrackParamEstimationAlg(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TracccDeviceTrackParamEstimationAlg", std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.detector == nullptr) {
    throw std::invalid_argument(
        "TracccDeviceTrackParamEstimationAlg: detector is not set");
  }
  if (m_cfg.magneticField == nullptr) {
    throw std::invalid_argument(
        "TracccDeviceTrackParamEstimationAlg: magneticField is not set");
  }
  m_algorithm = m_cfg.detector->backend().makeSeedParameterEstimation(
      m_cfg.estimation, this->logger().clone("TracccTrackParamEstimation"));

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSpacepoints.initialize(m_cfg.inputSpacepoints);
  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
}

ProcessCode TracccDeviceTrackParamEstimationAlg::execute(
    const AlgorithmContext& ctx) const {
  const TracccDetector& detector = *m_cfg.detector;

  TracccDeviceTrackParameters parameters = (**m_algorithm)(
      m_cfg.magneticField->deviceField(), m_inputMeasurements(ctx),
      m_inputSpacepoints(ctx), m_inputSeeds(ctx));
  detector.backend().synchronize();

  ACTS_DEBUG("Estimated " << m_algorithm->copy().get_size(parameters)
                          << " track parameters");

  m_outputTrackParameters(ctx, std::move(parameters));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
