// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccSeqAlg.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

// This is a separate file that maintains all the cuda headers
#include "ActsExamples/Traccc/TracccChain.hpp"

#include <stdexcept>

namespace ActsExamples {

TracccSeqAlgorithm::TracccSeqAlgorithm(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TracccSeqAlgorithm", std::move(logger)), m_cfg(cfg) {
  if (m_cfg.detectorFile.empty()) {
    throw std::invalid_argument("Detector geometry file is not configured");
  }
  if (m_cfg.conditionsFile.empty()) {
    throw std::invalid_argument("Detector conditions file is not configured");
  }
  if (m_cfg.digitizationFile.empty()) {
    throw std::invalid_argument("Detector digitization file is not configured");
  }
  if (m_cfg.bfieldFile.empty()) {
    throw std::invalid_argument("Magnetic field file is not configured");
  }
  if (m_cfg.dataDirectory.empty()) {
    throw std::invalid_argument("Data directory is not configured");
  }

  m_chain = std::make_shared<TracccChain>(
      m_cfg.detectorFile, m_cfg.digitizationFile, m_cfg.conditionsFile,
      m_cfg.materialFile, m_cfg.gridFile, m_cfg.bfieldFile);
}

ProcessCode TracccSeqAlgorithm::execute(const AlgorithmContext& ctx) const {
  const std::size_t eventId = ctx.eventNumber;

  ACTS_INFO("Processing event " << eventId);

  EventResult result = processEvent(m_chain, m_cfg.dataDirectory, eventId);

  ACTS_INFO("Event information:");
  ACTS_INFO("read " << result.n_cells << " cells,");
  ACTS_INFO("created " << result.n_measurements << " measurements,");
  ACTS_INFO("created " << result.n_spacepoints << "spacepoints,");
  ACTS_INFO("reconstructed " << result.n_seeds << " seeds,");
  ACTS_INFO("found " << result.n_found_tracks << " tracks,");
  ACTS_INFO("resolved " << result.n_resolved_tracks << " tracks,");
  ACTS_INFO("and fitted " << result.n_fitted_tracks << " tracks.");

  return ProcessCode::SUCCESS;
}

ProcessCode TracccSeqAlgorithm::finalize() {
  ACTS_INFO("Finalizing traccc GPU sequence algorithm");
  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
