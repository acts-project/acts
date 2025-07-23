// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/MuonSpacePointDigitizer.hpp"

#include "ActsExamples/EventData/MuonSpacePoint.hpp"

namespace ActsExamples {
MuonSpacePointDigitizer::MuonSpacePointDigitizer(const Config& cfg,
                                                 Acts::Logging::Level lvl)
    : IAlgorithm("MuonSpacePointDigitizer", lvl),
      m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("MuonSpacePointDigitizer", lvl)) {}

ProcessCode MuonSpacePointDigitizer::initialize() {
  return ProcessCode::SUCCESS;
}

ProcessCode MuonSpacePointDigitizer::execute(
    const AlgorithmContext& ctx) const {
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
