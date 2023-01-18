// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometryBuilder.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

#include <functional>

Acts::TrackingGeometryBuilder::TrackingGeometryBuilder(
    const Acts::TrackingGeometryBuilder::Config& cgbConfig,
    std::unique_ptr<const Logger> logger)
    : m_cfg(), m_logger(std::move(logger)) {
  setConfiguration(cgbConfig);
}

void Acts::TrackingGeometryBuilder::setConfiguration(
    const Acts::TrackingGeometryBuilder::Config& cgbConfig) {
  if (cgbConfig.trackingVolumeBuilders.empty()) {
    throw std::invalid_argument("Invalid configuration: no volume builders");
  }
  m_cfg = cgbConfig;
}

void Acts::TrackingGeometryBuilder::setLogger(
    std::unique_ptr<const Logger> newLogger) {
  m_logger = std::move(newLogger);
}

std::unique_ptr<const Acts::TrackingGeometry>
Acts::TrackingGeometryBuilder::trackingGeometry(
    const GeometryContext& gctx) const {
  MutableTrackingVolumePtr highestVolume = nullptr;
  // loop over the builders and wrap one around the other
  for (auto& volumeBuilder : m_cfg.trackingVolumeBuilders) {
    // assign a new highest volume (and potentially wrap around the given
    // highest volume so far)
    highestVolume = volumeBuilder(gctx, highestVolume, nullptr);
  }

  // create the TrackingGeometry & decorate it with the material
  if (highestVolume) {
    return std::make_unique<TrackingGeometry>(
        highestVolume,
        m_cfg.materialDecorator ? m_cfg.materialDecorator.get() : nullptr,
        *m_cfg.geometryIdentifierHook);
  } else {
    throw std::runtime_error(
        "Unable to construct tracking geometry: no tracking volume");
  }
}
