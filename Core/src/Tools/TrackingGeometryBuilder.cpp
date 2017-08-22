// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingGeometryBuilder.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Tools/TrackingGeometryBuilder.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"

Acts::TrackingGeometryBuilder::TrackingGeometryBuilder(
    const Acts::TrackingGeometryBuilder::Config& cgbConfig,
    std::unique_ptr<const Logger>                logger)
  : m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(cgbConfig);
}

void
Acts::TrackingGeometryBuilder::setConfiguration(
    const Acts::TrackingGeometryBuilder::Config& cgbConfig)
{
  // @todo check consistency
  // copy the configuration
  m_cfg = cgbConfig;
}

void
Acts::TrackingGeometryBuilder::setLogger(
    std::unique_ptr<const Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

std::unique_ptr<const Acts::TrackingGeometry>
Acts::TrackingGeometryBuilder::trackingGeometry() const
{
  // the return geometry with the highest volume
  std::unique_ptr<const TrackingGeometry> trackingGeometry;
  MutableTrackingVolumePtr                highestVolume = nullptr;
  // loop over the builders and wrap one around the other
  // -----------------------------
  for (auto& volumeBuilder : m_cfg.trackingVolumeBuilders) {
    // assign a new highest volume (and potentially wrap around the given
    // highest volume so far)
    highestVolume = volumeBuilder->trackingVolume(highestVolume);
  }  // --------------------------------------------------------------------------------

  // create the TrackingGeoemtry
  if (highestVolume)
    trackingGeometry.reset(new TrackingGeometry(highestVolume));
  // return the geometry to the service
  return (trackingGeometry);
}
