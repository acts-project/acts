// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DispatchAlgorithms/PatternDispatchAlgorithm.hpp"

namespace ActsExamples {

PatternDispatchAlgorithm::PatternDispatchAlgorithm(Config config, Acts::Logging::Level level)
    : IAlgorithm("PatternDispatchAlgorithm", level), m_cfg(std::move(config)) {}

ProcessCode PatternDispatchAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Retrieve the input data from the context
  const auto& measurementsIn = m_inputMeasurements(ctx);



  return ProcessCode::SUCCESS;
}

}