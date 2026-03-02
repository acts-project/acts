// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Svg/SvgTrackingGeometryWriter.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

namespace ActsExamples {

SvgTrackingGeometryWriter::SvgTrackingGeometryWriter(
    const SvgTrackingGeometryWriter::Config& config, Acts::Logging::Level level)
    : m_logger{Acts::getDefaultLogger(name(), level)}, m_cfg(config) {}

std::string SvgTrackingGeometryWriter::name() const {
  return "SvgTrackingGeometryWriter";
}

ProcessCode SvgTrackingGeometryWriter::write(
    const AlgorithmContext& context, const Acts::TrackingGeometry& tGeometry) {
  ACTS_DEBUG(">>Svg: Writer for TrackingGeometry object called.");

  m_writeMutex.lock();

  auto geometrySheets = ActsPlugins::Svg::TrackingGeometryConverter::convert(
      context.geoContext, tGeometry, m_cfg.converterOptions);

  // Write them out
  for (const auto& sheet : geometrySheets) {
    ActsPlugins::Svg::toFile({sheet},
                             joinPaths(m_cfg.outputDir, sheet._id + ".svg"));
  }
  // Successfully done
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
