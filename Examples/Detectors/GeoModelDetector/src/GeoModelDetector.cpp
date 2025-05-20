// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"

#include "GeoModelKernel/throwExcept.h"
namespace ActsExamples {

GeoModelDetector::GeoModelDetector(const Config& cfg)
    : Detector(Acts::getDefaultLogger("GeoModelDetector", cfg.logLevel)),
      m_cfg(cfg) {
  if (!m_cfg.geoModelTree.worldVolume) {
    m_cfg.geoModelTree = Acts::GeoModelReader::readFromDb(m_cfg.path);
  }
  if (!m_cfg.geoModelTree.worldVolume) {
    THROW_EXCEPTION("Failed to load geometry from '" << m_cfg.path << "'");
  }
  //// Dummy tracking geometry for the moments

  auto worldTrkVol = std::make_unique<Acts::TrackingVolume>(
      Acts::Transform3::Identity(),
      std::make_unique<Acts::CuboidVolumeBounds>(10, 10, 10), "SmallWorld");

  // m_trackingGeometry =
  // std::make_unique<Acts::TrackingGeometry>(std::move(worldTrkVol));

  // m_nominalGeometryContext = Acts::GeometryContext();

  // m_trackingGeometry =
  // m_cfg.trackingGeometryBuilder.trackingGeometry(m_nominalGeometryContext);
}

}  // namespace ActsExamples
