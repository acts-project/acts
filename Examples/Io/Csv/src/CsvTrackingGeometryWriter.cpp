// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvTrackingGeometryWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Plugins/Digitization/CartesianSegmentation.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

using namespace ActsExamples;

CsvTrackingGeometryWriter::CsvTrackingGeometryWriter(
    const CsvTrackingGeometryWriter::Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_world(nullptr),
      m_logger(Acts::getDefaultLogger("CsvTrackingGeometryWriter", level))

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

using SurfaceWriter = dfe::NamedTupleCsvWriter<SurfaceData>;
using SurfaceGridWriter = dfe::NamedTupleCsvWriter<SurfaceGridData>;
using LayerVolumeWriter = dfe::NamedTupleCsvWriter<LayerVolumeData>;
using BoundarySurface = Acts::BoundarySurfaceT<Acts::TrackingVolume>;

/// Write a single surface.
void fillSurfaceData(SurfaceData& data, const Acts::Surface& surface,
                     const Acts::GeometryContext& geoCtx) noexcept(false) {
  // encoded and partially decoded geometry identifier
  data.geometry_id = surface.geometryId().value();
  data.volume_id = surface.geometryId().volume();
  data.boundary_id = surface.geometryId().boundary();
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

  std::array<float*, 7> dataBoundParameters = {
      &data.bound_param0, &data.bound_param1, &data.bound_param2,
      &data.bound_param3, &data.bound_param4, &data.bound_param5,
      &data.bound_param6};

  const auto& bounds = surface.bounds();
  data.bounds_type = static_cast<int>(bounds.type());
  auto boundValues = bounds.values();

  if (boundValues.size() > dataBoundParameters.size()) {
    throw std::invalid_argument(
        "Bound types with too many parameters. Should never happen.");
  }

  for (size_t ipar = 0; ipar < boundValues.size(); ++ipar) {
    (*dataBoundParameters[ipar]) = boundValues[ipar];
  }
}

/// Write a single surface.
void writeSurface(SurfaceWriter& sfWriter, const Acts::Surface& surface,
                  const Acts::GeometryContext& geoCtx) {
  SurfaceData data;
  fillSurfaceData(data, surface, geoCtx);
  sfWriter.append(data);
}

/// Write a single surface.
void writeBoundarySurface(SurfaceWriter& writer,
                          const BoundarySurface& bsurface,
                          const Acts::GeometryContext& geoCtx) {
  SurfaceData data;
  fillSurfaceData(data, bsurface.surfaceRepresentation(), geoCtx);
  writer.append(data);
}

/// Write all child surfaces and descend into confined volumes.
void writeVolume(SurfaceWriter& sfWriter, SurfaceGridWriter& sfGridWriter,
                 LayerVolumeWriter& lvWriter,
                 const Acts::TrackingVolume& volume, bool writeSensitive,
                 bool writeBoundary, bool writeSurfaceGrid,
                 bool writeLayerVolume, const Acts::GeometryContext& geoCtx) {
  // process all layers that are directly stored within this volume
  if (volume.confinedLayers()) {
    const auto& vTransform = volume.transform();

    // Get the values of the volume boundaries
    std::vector<Acts::ActsScalar> volumeBoundValues =
        volume.volumeBounds().values();
    std::vector<Acts::ActsScalar> lastBoundValues;
    if (volume.volumeBounds().type() == Acts::VolumeBounds::eCylinder) {
      auto vTranslation = vTransform.translation();
      // change volume Bound values to r min , r max , z min, z max, phi min ,
      // phi max
      volumeBoundValues = {
          volumeBoundValues[0],
          volumeBoundValues[1],
          vTranslation.z() - volumeBoundValues[2],
          vTranslation.z() + volumeBoundValues[2],
          volumeBoundValues[4] - volumeBoundValues[3],
          volumeBoundValues[4] + volumeBoundValues[3],
      };
      lastBoundValues = volumeBoundValues;
    }

    /// Helper method for layer volume writing
    ///
    /// @param lv the layer volume to be written
    /// @param transform the layer transform
    /// @param representingBoundValues [in,out] the bound values
    /// @param last is the last layer
    auto writeCylinderLayerVolume =
        [&](const Acts::Layer& lv, const Acts::Transform3& transform,
            std::vector<Acts::ActsScalar>& representingBoundValues,
            bool last) -> void {
      // The layer volume to be written
      LayerVolumeData lvDims;
      lvDims.geometry_id = lv.geometryId().value();
      lvDims.volume_id = lv.geometryId().volume();
      lvDims.layer_id = lv.geometryId().layer();
      bool isCylinderLayer = (lv.surfaceRepresentation().bounds().type() ==
                              Acts::SurfaceBounds::eCylinder);

      auto lTranslation = transform.translation();
      // Change volume Bound values to r min , r max , z min, z max, phi min ,
      // phi max
      representingBoundValues = {
          representingBoundValues[0],
          representingBoundValues[1],
          lTranslation.z() - representingBoundValues[2],
          lTranslation.z() + representingBoundValues[2],
          representingBoundValues[4] - representingBoundValues[3],
          representingBoundValues[4] + representingBoundValues[3]};

      // Synchronize
      lvDims.min_v0 =
          isCylinderLayer ? representingBoundValues[0] : volumeBoundValues[0];
      lvDims.max_v0 =
          isCylinderLayer ? representingBoundValues[1] : volumeBoundValues[1];
      lvDims.min_v1 =
          isCylinderLayer ? volumeBoundValues[2] : representingBoundValues[2];
      lvDims.max_v1 =
          isCylinderLayer ? volumeBoundValues[3] : representingBoundValues[3];
      lvDims.min_v2 = representingBoundValues[4];
      lvDims.max_v2 = representingBoundValues[5];

      // Write the prior navigation layer
      LayerVolumeData nlvDims;
      nlvDims.geometry_id = lv.geometryId().value();
      nlvDims.volume_id = lv.geometryId().volume();
      nlvDims.layer_id = lv.geometryId().layer() - 1;
      if (isCylinderLayer) {
        nlvDims.min_v0 = lastBoundValues[0];
        nlvDims.max_v0 = representingBoundValues[0];
        nlvDims.min_v1 = lastBoundValues[2];
        nlvDims.max_v1 = lastBoundValues[3];
        // Reset the r min boundary
        lastBoundValues[0] = representingBoundValues[1];
      } else {
        nlvDims.min_v0 = lastBoundValues[0];
        nlvDims.max_v0 = lastBoundValues[1];
        nlvDims.min_v1 = lastBoundValues[2];
        nlvDims.max_v1 = representingBoundValues[2];
        // Reset the r min boundary
        lastBoundValues[2] = representingBoundValues[3];
      }
      nlvDims.min_v2 = representingBoundValues[4];
      nlvDims.max_v2 = representingBoundValues[5];
      lvWriter.append(nlvDims);
      // Write the volume dimensions for the sensitive layer
      lvWriter.append(lvDims);

      // Write last if needed
      if (last) {
        // Write the last navigation layer volume
        LayerVolumeData llvDims;
        llvDims.geometry_id = lv.geometryId().value();
        llvDims.volume_id = lv.geometryId().volume();
        llvDims.layer_id = lv.geometryId().layer() + 1;
        llvDims.min_v0 = lastBoundValues[0];
        llvDims.max_v0 = lastBoundValues[1];
        llvDims.min_v1 = lastBoundValues[2];
        llvDims.max_v1 = lastBoundValues[3];
        llvDims.min_v2 = representingBoundValues[4];
        llvDims.max_v2 = representingBoundValues[5];
        // Close up volume
        lvWriter.append(llvDims);
      }
    };

    unsigned int layerIdx = 0;
    const auto& layers = volume.confinedLayers()->arrayObjects();
    for (auto layer : layers) {
      // We jump navigation layers - they will be written last
      if (layer->layerType() == Acts::navigation) {
        ++layerIdx;
        // For a correct layer volume setup, we need the navigation layers
        if (writeLayerVolume) {
          writeSurface(sfWriter, layer->surfaceRepresentation(), geoCtx);
        }
        continue;

      } else if (layer->surfaceArray()) {
        // Get the representing volume
        const auto* rVolume = layer->representingVolume();
        // Write the layer volume
        if (rVolume != nullptr and writeLayerVolume) {
          // Get the values of the representing volume
          std::vector<Acts::ActsScalar> representingBoundValues =
              rVolume->volumeBounds().values();
          if (rVolume->volumeBounds().type() == Acts::VolumeBounds::eCylinder) {
            bool last = (layerIdx + 2 ==
                         volume.confinedLayers()->arrayObjects().size());
            writeCylinderLayerVolume(*layer, rVolume->transform(),
                                     representingBoundValues, last);
          }
        }

        auto sfArray = layer->surfaceArray();
        // Write the surface grid
        if (writeSurfaceGrid) {
          SurfaceGridData sfGrid;
          sfGrid.geometry_id = layer->geometryId().value();
          sfGrid.volume_id = layer->geometryId().volume();
          sfGrid.layer_id = layer->geometryId().layer();

          // Draw the grid itself
          auto binning = sfArray->binningValues();
          auto axes = sfArray->getAxes();
          if (not binning.empty() and binning.size() == 2 and
              axes.size() == 2) {
            auto loc0Values = axes[0]->getBinEdges();
            sfGrid.nbins_loc0 = loc0Values.size();
            sfGrid.min_loc0 = loc0Values[0];
            sfGrid.max_loc0 = loc0Values[loc0Values.size() - 1];
            auto loc1Values = axes[1]->getBinEdges();
            sfGrid.nbins_loc1 = loc1Values.size();
            sfGrid.min_loc1 = loc1Values[0];
            sfGrid.max_loc1 = loc1Values[loc1Values.size() - 1];
          }

          sfGridWriter.append(sfGrid);
        }

        if (writeSensitive) {
          for (auto surface : sfArray->surfaces()) {
            if (surface) {
              writeSurface(sfWriter, *surface, geoCtx);
            }
          }
        }
      } else {
        writeSurface(sfWriter, layer->surfaceRepresentation(), geoCtx);
        // Layer volume definition of passive layers
        if (writeLayerVolume) {
          if (volume.volumeBounds().type() == Acts::VolumeBounds::eCylinder) {
            // Only one passive layer present - write volume directly
            if (layers.size() == 3) {
              LayerVolumeData plvDims;
              plvDims.geometry_id = layer->geometryId().value();
              plvDims.volume_id = layer->geometryId().volume();
              plvDims.layer_id = layer->geometryId().layer();
              plvDims.min_v0 = volumeBoundValues[0];
              plvDims.max_v0 = volumeBoundValues[1];
              plvDims.min_v1 = volumeBoundValues[2];
              plvDims.max_v1 = volumeBoundValues[3];
              lvWriter.append(plvDims);
            }
          }
        }
      }
      ++layerIdx;
    }  // enf of layer loop

    // This is a navigation volume, write the boundaries
    if (writeBoundary) {
      for (auto bsurface : volume.boundarySurfaces()) {
        writeBoundarySurface(sfWriter, *bsurface, geoCtx);
      }
    }
  }
  // step down into hierarchy to process all child volumnes
  if (volume.confinedVolumes()) {
    for (auto confined : volume.confinedVolumes()->arrayObjects()) {
      writeVolume(sfWriter, sfGridWriter, lvWriter, *confined.get(),
                  writeSensitive, writeBoundary, writeSurfaceGrid,
                  writeLayerVolume, geoCtx);
    }
  }
}
}  // namespace

ProcessCode CsvTrackingGeometryWriter::write(const AlgorithmContext& ctx) {
  if (not m_cfg.writePerEvent) {
    return ProcessCode::SUCCESS;
  }

  SurfaceWriter sfWriter(
      perEventFilepath(m_cfg.outputDir, "detectors.csv", ctx.eventNumber),
      m_cfg.outputPrecision);

  SurfaceGridWriter sfGridWriter(
      perEventFilepath(m_cfg.outputDir, "surface-grids.csv", ctx.eventNumber),
      m_cfg.outputPrecision);

  LayerVolumeWriter lvWriter(
      perEventFilepath(m_cfg.outputDir, "layer-volumes.csv", ctx.eventNumber),
      m_cfg.outputPrecision);

  writeVolume(sfWriter, sfGridWriter, lvWriter, *m_world, m_cfg.writeSensitive,
              m_cfg.writeBoundary, m_cfg.writeSurfaceGrid,
              m_cfg.writeLayerVolume, ctx.geoContext);
  return ProcessCode::SUCCESS;
}

ProcessCode CsvTrackingGeometryWriter::endRun() {
  SurfaceWriter sfWriter(joinPaths(m_cfg.outputDir, "detectors.csv"),
                         m_cfg.outputPrecision);
  SurfaceGridWriter sfGridWriter(
      joinPaths(m_cfg.outputDir, "surface-grids.csv"), m_cfg.outputPrecision);

  LayerVolumeWriter lvWriter(joinPaths(m_cfg.outputDir, "layer-volumes.csv"),
                             m_cfg.outputPrecision);

  writeVolume(sfWriter, sfGridWriter, lvWriter, *m_world, m_cfg.writeSensitive,
              m_cfg.writeBoundary, m_cfg.writeSurfaceGrid,
              m_cfg.writeLayerVolume, Acts::GeometryContext());
  return ProcessCode::SUCCESS;
}
