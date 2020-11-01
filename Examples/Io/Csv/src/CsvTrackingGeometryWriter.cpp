// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvTrackingGeometryWriter.hpp"

#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Plugins/Digitization/CartesianSegmentation.hpp>
#include <Acts/Plugins/Digitization/DigitizationModule.hpp>
#include <Acts/Plugins/Identification/IdentifiedDetectorElement.hpp>
#include <Acts/Surfaces/DiscBounds.hpp>
#include <Acts/Surfaces/PlanarBounds.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Units.hpp>

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvDataFormats.hpp"

using namespace ActsExamples;

CsvTrackingGeometryWriter::CsvTrackingGeometryWriter(
    const CsvTrackingGeometryWriter::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_world(nullptr),
      m_logger(Acts::getDefaultLogger("CsvTrackingGeometryWriter", lvl))

{
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  m_world = m_cfg.trackingGeometry->highestTrackingVolume();
  if (not m_world) {
    throw std::invalid_argument("Could not identify the world volume");
  }
}

std::string CsvTrackingGeometryWriter::name() const {
  return "CsvTrackingGeometryWriter";
}

namespace {

using SegmentWriter = dfe::NamedTupleCsvWriter<SegmentData>;

using SurfaceWriter = dfe::NamedTupleCsvWriter<SurfaceData>;

/// Write a single surface.
void writeSurface(SegmentWriter& swriter, SurfaceWriter& writer,
                  const Acts::Surface& surface,
                  const Acts::GeometryContext& geoCtx) {
  SurfaceData data;

  // encoded and partially decoded geometry identifier
  data.geometry_id = surface.geometryId().value();
  data.volume_id = surface.geometryId().volume();
  data.layer_id = surface.geometryId().layer();
  data.module_id = surface.geometryId().sensitive();
  // center position
  auto center = surface.center(geoCtx);
  data.cx = center.x() / Acts::UnitConstants::mm;
  data.cy = center.y() / Acts::UnitConstants::mm;
  data.cz = center.z() / Acts::UnitConstants::mm;
  // rotation matrix components are unit-less
  auto transform = surface.transform(geoCtx);
  data.rot_xu = transform(0, 0);
  data.rot_xv = transform(0, 1);
  data.rot_xw = transform(0, 2);
  data.rot_yu = transform(1, 0);
  data.rot_yv = transform(1, 1);
  data.rot_yw = transform(1, 2);
  data.rot_zu = transform(2, 0);
  data.rot_zv = transform(2, 1);
  data.rot_zw = transform(2, 2);

  // module thickness
  if (surface.associatedDetectorElement()) {
    const auto* detElement =
        dynamic_cast<const Acts::IdentifiedDetectorElement*>(
            surface.associatedDetectorElement());
    if (detElement) {
      data.module_t = detElement->thickness() / Acts::UnitConstants::mm;
    }
  }

  // bounds and pitch (if available)
  std::vector<Acts::Vector2D> vertices;
  const Acts::SurfaceBounds& bounds = surface.bounds();
  const auto* planarBounds = dynamic_cast<const Acts::PlanarBounds*>(&bounds);
  if (planarBounds) {
    // Get the vertices
    vertices = planarBounds->vertices(1);

    // extract limits from value store
    auto boundValues = surface.bounds().values();
    if (boundValues.size() == 2) {
      data.module_minhu = boundValues[0] / Acts::UnitConstants::mm;
      data.module_minhu = boundValues[0] / Acts::UnitConstants::mm;
      data.module_minhu = boundValues[1] / Acts::UnitConstants::mm;
    } else if (boundValues.size() == 3) {
      data.module_minhu = boundValues[0] / Acts::UnitConstants::mm;
      data.module_minhu = boundValues[0] / Acts::UnitConstants::mm;
      data.module_minhu = boundValues[1] / Acts::UnitConstants::mm;
    }
    // get the pitch from the digitization module
    const auto* detElement =
        dynamic_cast<const Acts::IdentifiedDetectorElement*>(
            surface.associatedDetectorElement());
    if (detElement and detElement->digitizationModule()) {
      auto dModule = detElement->digitizationModule();
      // dynamic_cast to CartesianSegmentation
      const auto* cSegmentation =
          dynamic_cast<const Acts::CartesianSegmentation*>(
              &(dModule->segmentation()));
      if (cSegmentation) {
        auto pitch = cSegmentation->pitch();
        data.pitch_u = pitch.first / Acts::UnitConstants::mm;
        data.pitch_u = pitch.second / Acts::UnitConstants::mm;
      }
    }
  } else {
    const auto* discBounds = dynamic_cast<const Acts::DiscBounds*>(&bounds);
    // Get the vertices
    vertices = discBounds->vertices(72);
  }

  writer.append(data);

  if (not vertices.empty()) {
    // Close the loop
    vertices.push_back(vertices[0]);
    for (size_t iv = 1; iv < vertices.size(); ++iv) {
      Acts::Vector2D lstart = vertices[iv - 1];
      Acts::Vector2D lend = vertices[iv];
      Acts::Vector3D start =
          transform * Acts::Vector3D(lstart.x(), lstart.y(), 0.);
      Acts::Vector3D end = transform * Acts::Vector3D(lend.x(), lend.y(), 0.);
      SegmentData sdata;
      sdata.geometry_id = data.geometry_id;
      sdata.volume_id = data.volume_id;
      sdata.layer_id = data.layer_id;
      sdata.module_id = data.module_id;
      sdata.type = 'l';
      sdata.object = 's';
      sdata.p0 = start.x();
      sdata.p1 = start.y();
      sdata.p2 = end.x();
      sdata.p3 = end.y();
      swriter.append(sdata);
    }
  }
}

/// Write all child surfaces and descend into confined volumes.
void writeVolume(SegmentWriter& swriter, SurfaceWriter& writer,
                 const Acts::TrackingVolume& volume,
                 const Acts::GeometryContext& geoCtx) {
  // process all layers that are directly stored within this volume
  if (volume.confinedLayers()) {
    for (auto layer : volume.confinedLayers()->arrayObjects()) {
      // we jump navigation layers
      if (layer->layerType() == Acts::navigation) {
        continue;
      }
      // check for sensitive surfaces
      if (layer->surfaceArray()) {
        for (auto surface : layer->surfaceArray()->surfaces()) {
          if (surface) {
            writeSurface(swriter, writer, *surface, geoCtx);
          }
        }
      }
    }
  }
  // step down into hierarchy to process all child volumnes
  if (volume.confinedVolumes()) {
    for (auto confined : volume.confinedVolumes()->arrayObjects()) {
      writeVolume(swriter, writer, *confined.get(), geoCtx);
    }
  }
}
}  // namespace

ProcessCode CsvTrackingGeometryWriter::write(const AlgorithmContext& ctx) {
  if (not m_cfg.writePerEvent) {
    return ProcessCode::SUCCESS;
  }
  SegmentWriter swriter(joinPaths(m_cfg.outputDir, "segments.csv"),
                        m_cfg.outputPrecision);

  SurfaceWriter writer(joinPaths(m_cfg.outputDir, "detectors.csv"),
                       m_cfg.outputPrecision);
  writeVolume(swriter, writer, *m_world, Acts::GeometryContext());
  return ProcessCode::SUCCESS;
}

ProcessCode CsvTrackingGeometryWriter::endRun() {
  SegmentWriter swriter(joinPaths(m_cfg.outputDir, "segments.csv"),
                        m_cfg.outputPrecision);

  SurfaceWriter writer(joinPaths(m_cfg.outputDir, "detectors.csv"),
                       m_cfg.outputPrecision);
  writeVolume(swriter, writer, *m_world, Acts::GeometryContext());
  return ProcessCode::SUCCESS;
}
