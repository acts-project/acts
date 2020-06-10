// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename seedfinder_t>
SeedingAlgorithm<seedfinder_t>::SeedingAlgorithm(
    SeedingAlgorithm<seedfinder_t>::Config cfg, Acts::Logging::Level lvl)
    : BareAlgorithm(lvl), m_cfg(cfg) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing input space point collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing output seeds collection");
  }
}

template <typename seedfinder_t>
SeedingAlgorithm<seedfinder_t>::execute(const AlgorithmContext& ctx) const {}
