// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Plugins/Svg/SvgTrackingGeometryWriter.hpp"

#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>

#include <iostream>

ActsExamples::SvgTrackingGeometryWriter::SvgTrackingGeometryWriter(
    const ActsExamples::SvgTrackingGeometryWriter::Config& config,
    Acts::Logging::Level level)
    : m_logger{Acts::getDefaultLogger(name(), level)}, m_cfg(config) {}

std::string ActsExamples::SvgTrackingGeometryWriter::name() const {
  return "SvgTrackingGeometryWriter";
}

ActsExamples::ProcessCode ActsExamples::SvgTrackingGeometryWriter::write(
    const AlgorithmContext& context, const Acts::TrackingGeometry& tGeometry) {
  ACTS_DEBUG(">>Svg: Writer for TrackingGeometry Svgect called.");

  auto world = tGeometry.highestTrackingVolume();
  if (world) {
    write(context, *world);
  }
  return ActsExamples::ProcessCode::SUCCESS;
}

void ActsExamples::SvgTrackingGeometryWriter::write(
    const AlgorithmContext& context, const Acts::TrackingVolume& tVolume) {
  ACTS_DEBUG(">>Svg: Writer for TrackingVolume Svgect called.");

}
