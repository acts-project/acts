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
  // @todo check consistency
  // copy the configuration
  m_cfg = cgbConfig;
}

void Acts::TrackingGeometryBuilder::setLogger(
    std::unique_ptr<const Logger> newLogger) {
  m_logger = std::move(newLogger);
}

std::unique_ptr<const Acts::TrackingGeometry>
Acts::TrackingGeometryBuilder::trackingGeometry(
    const GeometryContext& gctx) const {
  // the return geometry with the highest volume
  std::unique_ptr<const TrackingGeometry> trackingGeometry;
  MutableTrackingVolumePtr highestVolume = nullptr;
  // loop over the builders and wrap one around the other
  // -----------------------------
  for (auto& volumeBuilder : m_cfg.trackingVolumeBuilders) {
    // assign a new highest volume (and potentially wrap around the given
    // highest volume so far)
    highestVolume = volumeBuilder(gctx, highestVolume, nullptr);
  }  // --------------------------------------------------------------------------------

  // create the TrackingGeometry & decorate it with the material
  if (highestVolume) {
    // first check if we have material to get
    const IMaterialDecorator* materialDecorator =
        m_cfg.materialDecorator ? m_cfg.materialDecorator.get() : nullptr;
    // build and set the TrackingGeometry
    trackingGeometry.reset(new TrackingGeometry(
        highestVolume, materialDecorator, m_cfg.geometryIdentifierHook));
  }
  // return the geometry to the service
  return (trackingGeometry);
}
