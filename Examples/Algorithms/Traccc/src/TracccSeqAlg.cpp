// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccSeqAlg.hpp"

namespace ActsExamples::Traccc {

TracccSeqAlg::TracccSeqAlg(Config cfg, Acts::Logging::Level logLevel)
    : IAlgorithm("TracccSeqAlg", logLevel),
      m_cfg(std::move(cfg)),
      m_chain(m_cfg.chain, logLevel) {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSpacepoints.initialize(m_cfg.inputSpacepoints);
}

ProcessCode TracccSeqAlg::execute(const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& spacepoints = m_inputSpacepoints(ctx);

  ACTS_DEBUG("Running traccc chain on event "
             << ctx.eventNumber << " with " << measurements.size()
             << " measurements and " << spacepoints.size() << " spacepoints.");

  EventResult result = m_chain(
      const_cast<traccc::edm::measurement_collection::host&>(measurements),
      const_cast<traccc::edm::spacepoint_collection::host&>(spacepoints));

  ACTS_INFO("Event information:");

  ACTS_INFO("received " << result.n_measurements << " measurements,");
  ACTS_INFO("received " << result.n_spacepoints << " spacepoints,");
  ACTS_INFO("reconstructed " << result.n_seeds << " seeds,");
  ACTS_INFO("found " << result.n_track_candidates << " tracks,");
  ACTS_INFO("and fitted " << result.n_fitted_tracks << " tracks.");

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples::Traccc
