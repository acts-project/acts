// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Svg/SvgTrackingGeometryWriter.hpp"

#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Plugins/ActSVG/LayerSvgConverter.hpp>

#include <iostream>
#include <string>

ActsExamples::SvgTrackingGeometryWriter::SvgTrackingGeometryWriter(
    const ActsExamples::SvgTrackingGeometryWriter::Config& config,
    Acts::Logging::Level level)
    : m_logger{Acts::getDefaultLogger(name(), level)}, m_cfg(config) {}

std::string ActsExamples::SvgTrackingGeometryWriter::name() const {
  return "SvgTrackingGeometryWriter";
}

ActsExamples::ProcessCode ActsExamples::SvgTrackingGeometryWriter::write(
    const AlgorithmContext& context, const Acts::TrackingGeometry& tGeometry) {
  ACTS_DEBUG(">>Svg: Writer for TrackingGeometry object called.");

  m_writeMutex.lock();

  auto world = tGeometry.highestTrackingVolume();
  if (world) {
    write(context, *world);
  }

  Acts::Svg::toFile(m_xyCrossection, m_cfg.baseName + "crosssection_xy.svg");
  Acts::Svg::toFile(m_zrCrossection, m_cfg.baseName + "crosssection_zr.svg");

  return ActsExamples::ProcessCode::SUCCESS;
}

void ActsExamples::SvgTrackingGeometryWriter::write(
    const AlgorithmContext& context, const Acts::TrackingVolume& tVolume) {
  if (tVolume.confinedLayers() != nullptr) {
    for (const auto& layer : tVolume.confinedLayers()->arrayObjects()) {
      if (layer->surfaceArray() != nullptr) {
        Acts::GeometryIdentifier geoID = layer->geometryId();
        std::string layerName = m_cfg.baseName + "vol_" +
                                std::to_string(geoID.volume()) + "_layer_" +
                                std::to_string(geoID.layer());

        // Get the layer sheets
        auto layerSheet = Acts::Svg::layerSheets(
            context.geoContext, *layer, layerName, m_cfg.svgStyle,
            {m_cfg.zViewRangeXY, m_cfg.phiViewRangeRZ});
        // Write the individual sheets
        for (unsigned int is = 0; is < layerSheet.size(); ++is) {
          Acts::Svg::toFile({layerSheet[is]}, layerSheet[is]._id + ".svg");
        }
        // Collect the xy views
        if (layerSheet[2].is_defined() and
            layer->surfaceRepresentation().type() == Acts::Surface::Cylinder) {
          m_xyCrossection.push_back(layerSheet[2]);
        }
        // Collect the zr views
        if (layerSheet[3].is_defined()) {
          m_zrCrossection.push_back(layerSheet[3]);
        }
      }
    }
  }

  if (tVolume.confinedVolumes() != nullptr) {
    for (const auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      write(context, *sVolume);
    }
  }
}
