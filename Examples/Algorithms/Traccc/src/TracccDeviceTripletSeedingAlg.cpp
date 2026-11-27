// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccDeviceTripletSeedingAlg.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>

namespace ActsExamples {

TracccDeviceTripletSeedingAlg::TracccDeviceTripletSeedingAlg(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TracccDeviceTripletSeedingAlg", std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.detector == nullptr) {
    throw std::invalid_argument(
        "TracccDeviceTripletSeedingAlg: detector is not set");
  }

  m_cfg.seedFinder.setup();
  const traccc::spacepoint_grid_config gridConfig(m_cfg.seedFinder);

  m_algorithm = m_cfg.detector->backend().makeTripletSeeding(
      m_cfg.seedFinder, gridConfig, m_cfg.seedFilter,
      this->logger().clone("TracccTripletSeeding"));

  m_inputSpacepoints.initialize(m_cfg.inputSpacepoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
}

ProcessCode TracccDeviceTripletSeedingAlg::execute(
    const AlgorithmContext& ctx) const {
  const TracccDeviceSpacepoints& spacepoints = m_inputSpacepoints(ctx);

  TracccDeviceSeeds seeds = (**m_algorithm)(spacepoints);
  m_cfg.detector->backend().synchronize();

  ACTS_DEBUG("Found " << m_algorithm->copy().get_size(seeds) << " seeds");

  m_outputSeeds(ctx, std::move(seeds));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
