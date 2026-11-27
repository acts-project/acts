// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccDeviceSpacepointFormationAlg.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>

namespace ActsExamples {

TracccDeviceSpacepointFormationAlg::TracccDeviceSpacepointFormationAlg(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TracccDeviceSpacepointFormationAlg", std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.detector == nullptr) {
    throw std::invalid_argument(
        "TracccDeviceSpacepointFormationAlg: detector is not set");
  }
  m_algorithm = m_cfg.detector->backend().makeSpacepointFormation(
      this->logger().clone("TracccSpacepointFormation"));

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputSpacepoints.initialize(m_cfg.outputSpacepoints);
}

ProcessCode TracccDeviceSpacepointFormationAlg::execute(
    const AlgorithmContext& ctx) const {
  const TracccDeviceMeasurements& measurements = m_inputMeasurements(ctx);
  const TracccDetector& detector = *m_cfg.detector;

  TracccDeviceSpacepoints spacepoints =
      (**m_algorithm)(detector.deviceDetector(), measurements);
  detector.backend().synchronize();

  ACTS_DEBUG("Formed " << m_algorithm->copy().get_size(spacepoints)
                       << " spacepoints");

  m_outputSpacepoints(ctx, std::move(spacepoints));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
