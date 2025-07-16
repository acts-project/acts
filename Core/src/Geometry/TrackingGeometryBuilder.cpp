// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometryBuilder.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"

#include <functional>
#include <stdexcept>
#include <utility>

namespace Acts {

TrackingGeometryBuilder::TrackingGeometryBuilder(
    const TrackingGeometryBuilder::Config& cgbConfig,
    std::unique_ptr<const Logger> logger)
    : m_cfg(), m_logger(std::move(logger)) {
  setConfiguration(cgbConfig);
}

const TrackingGeometryBuilder::Config&
TrackingGeometryBuilder::getConfiguration() const {
  return m_cfg;
}

void TrackingGeometryBuilder::setConfiguration(
    const TrackingGeometryBuilder::Config& cgbConfig) {
  if (cgbConfig.trackingVolumeBuilders.empty()) {
    throw std::invalid_argument("Invalid configuration: no volume builders");
  }
  m_cfg = cgbConfig;
}

void TrackingGeometryBuilder::setLogger(
    std::unique_ptr<const Logger> newLogger) {
  m_logger = std::move(newLogger);
}

std::unique_ptr<const TrackingGeometry>
TrackingGeometryBuilder::trackingGeometry(const GeometryContext& gctx) const {
  ACTS_DEBUG("Building tracking geometry");
  MutableTrackingVolumePtr highestVolume = nullptr;
  // loop over the builders and wrap one around the other
  for (auto& volumeBuilder : m_cfg.trackingVolumeBuilders) {
    // assign a new highest volume (and potentially wrap around the given
    // highest volume so far)
    auto volume = volumeBuilder(gctx, highestVolume, nullptr);
    if (!volume) {
      ACTS_INFO(
          "Received nullptr volume from builder, keeping previous highest "
          "volume");
    } else {
      highestVolume = std::move(volume);
    }
  }

  // create the TrackingGeometry & decorate it with the material
  if (highestVolume) {
    return std::make_unique<TrackingGeometry>(
        highestVolume,
        m_cfg.materialDecorator ? m_cfg.materialDecorator.get() : nullptr,
        *m_cfg.geometryIdentifierHook, logger());
  } else {
    throw std::runtime_error(
        "Unable to construct tracking geometry: no tracking volume");
  }
}

}  // namespace Acts
