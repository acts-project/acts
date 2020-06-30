// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Obj/ObjTrackingGeometryWriter.hpp"

#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Visualization/GeometryView.hpp>
#include <Acts/Visualization/ObjVisualization.hpp>
#include <iostream>

FW::Obj::ObjTrackingGeometryWriter::ObjTrackingGeometryWriter(
    const FW::Obj::ObjTrackingGeometryWriter::Config& cfg)
    : m_cfg(cfg) {}

std::string FW::Obj::ObjTrackingGeometryWriter::name() const {
  return m_cfg.name;
}

FW::ProcessCode FW::Obj::ObjTrackingGeometryWriter::write(
    const AlgorithmContext& context, const Acts::TrackingGeometry& tGeometry) {
  ACTS_DEBUG(">>Obj: Writer for TrackingGeometry object called.");

  auto world = tGeometry.highestTrackingVolume();
  if (world) {
    write(context, *world);
  }
  return FW::ProcessCode::SUCCESS;
}

void FW::Obj::ObjTrackingGeometryWriter::write(
    const AlgorithmContext& context, const Acts::TrackingVolume& tVolume) {
  ACTS_DEBUG(">>Obj: Writer for TrackingVolume object called.");

  Acts::ObjVisualization objVis(m_cfg.outputPrecision, m_cfg.outputScalor);

  Acts::GeometryView::drawTrackingVolume(
      objVis, tVolume, context.geoContext, m_cfg.containerView,
      m_cfg.volumeView, m_cfg.sensitiveView, m_cfg.passiveView, m_cfg.gridView);
}
