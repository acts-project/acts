// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/BlueprintDetectorBuilder.hpp"

Acts::BlueprintDetectorBuilder::BlueprintDetectorBuilder(
    const Acts::BlueprintDetectorBuilder::Config& cfg,
    std::unique_ptr<const Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  m_cfg.protoDetector.harmonize(false);
}

std::shared_ptr<Acts::Experimental::Detector>
Acts::BlueprintDetectorBuilder::construct(const GeometryContext& gctx) const {
  return nullptr;
}
