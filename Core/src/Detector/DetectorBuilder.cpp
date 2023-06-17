// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/DetectorBuilder.hpp"

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"

#include <stdexcept>

Acts::Experimental::DetectorBuilder::DetectorBuilder(
    const Acts::Experimental::DetectorBuilder::Config& cfg,
    std::unique_ptr<const Acts::Logger> mlogger)
    : IDetectorBuilder(), m_cfg(cfg), m_logger(std::move(mlogger)) {
  if (m_cfg.builder == nullptr) {
    throw std::invalid_argument(
        "DetectorBuilder: no top level builder defined.");
  }
}

std::shared_ptr<const Acts::Experimental::Detector>
Acts::Experimental::DetectorBuilder::construct(
    const GeometryContext& gctx) const {
  // Screen printout of the auxilliary information
  if (not m_cfg.auxilliary.empty()) {
    ACTS_DEBUG(m_cfg.auxilliary);
  }
  ACTS_DEBUG("Building a detector with name " << m_cfg.name);

  auto [volumes, portals, roots] = m_cfg.builder->construct(gctx);

  return Detector::makeShared(m_cfg.name, std::move(roots.volumes),
                              std::move(roots.volumeFinder));
}
