// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DetectorBuilder.hpp"

#include "Acts/Geometry/detail/DetectorVolumeFinders.hpp"

Acts::Experimental::DetectorBuilder::DetectorBuilder(
    const Acts::Experimental::DetectorBuilder::Config& cfg,
    std::unique_ptr<const Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::shared_ptr<Acts::Experimental::Detector>
Acts::Experimental::DetectorBuilder::construct(
    const GeometryContext& gctx, const ProtoDetector& protoDetector) const {
  ACTS_DEBUG("Building Acts::Detector geometry from ProtoDetector '"
             << protoDetector.name << "'.");
  // Make a copy of the proto detector description, as we harmonize
  auto pDetector = protoDetector;
  pDetector.harmonize(false);
  // Get the world volume
  auto& worldVolume = pDetector.worldVolume;
  // Recursive blue print building
  DetectorBlock dBlock;
  worldVolume.blockBuilder(dBlock, gctx, m_cfg.logLevel);
  // Get the volumes that build this detector
  auto volumes = std::get<DetectorVolumes>(dBlock);
  // @TODO add differen volume finder generators

  // Return the detector
  return Detector::makeShared(pDetector.name, volumes, detail::tryAllVolumes());
}
