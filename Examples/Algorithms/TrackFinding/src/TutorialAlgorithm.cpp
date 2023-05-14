// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TutorialAlgorithm.hpp"

ActsExamples::TutorialAlgorithm::TutorialAlgorithm(const Config& cfg,
                                                   Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("TutorialAlgorithm", lvl), m_cfg(cfg) {}

ActsExamples::ProcessCode ActsExamples::TutorialAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  ACTS_INFO(m_cfg.message);

  return ActsExamples::ProcessCode::SUCCESS;
}
