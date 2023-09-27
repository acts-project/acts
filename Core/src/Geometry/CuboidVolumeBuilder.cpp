// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CuboidVolumeBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/detail/DefaultDetectorElementBase.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <limits>

std::shared_ptr<const Acts::Surface> Acts::CuboidVolumeBuilder::buildSurface(
    const GeometryContext& /*gctx*/,
    const CuboidVolumeBuilder::SurfaceConfig& cfg) const {
  std::shared_ptr<PlaneSurface> surface;

  // Build transformation
  Transform3 trafo(Transform3::Identity() * cfg.rotation);
  trafo.translation() = cfg.position;

  // Create and store surface
  if (cfg.detElementConstructor) {
    surface = Surface::makeShared<PlaneSurface>(
        cfg.rBounds,
        *(cfg.detElementConstructor(trafo, cfg.rBounds, cfg.thickness)));
  } else {
    surface = Surface::makeShared<PlaneSurface>(trafo, cfg.rBounds);
  }
  surface->assignSurfaceMaterial(cfg.surMat);
  return surface;
}

std::shared_ptr<const Acts::Layer> Acts::CuboidVolumeBuilder::buildLayer(
    const GeometryContext& gctx,
    Acts::CuboidVolumeBuilder::LayerConfig& cfg) const {
  if (cfg.surfaces.empty() && cfg.surfaceCfg.empty()) {
    throw std::runtime_error{
        "Neither surfaces nor config to build surfaces was provided. Cannot "
        "proceed"};
  }

  // Build the surface
  if (cfg.surfaces.empty()) {
    for (const auto& sCfg : cfg.surfaceCfg) {
      cfg.surfaces.push_back(buildSurface(gctx, sCfg));
    }
  }
  // Build transformation centered at the surface position
  Vector3 centroid{0., 0., 0.};

  for (const auto& surface : cfg.surfaces) {
    centroid += surface->transform(gctx).translation();
  }

  centroid /= cfg.surfaces.size();

  // In the case the layer configuration doesn't define the rotation of the
  // layer use the orientation of the first surface to define the layer rotation
  // in space.
  Transform3 trafo = Transform3::Identity();
  trafo.translation() = centroid;
  if (cfg.rotation) {
    trafo.linear() = *cfg.rotation;
  } else {
    trafo.linear() = cfg.surfaces.front()->transform(gctx).rotation();
  }

  LayerCreator::Config lCfg;
  lCfg.surfaceArrayCreator = std::make_shared<const SurfaceArrayCreator>();
  LayerCreator layerCreator(lCfg);
  ProtoLayer pl{gctx, cfg.surfaces};
  pl.envelope[binX] = cfg.envelopeX;
  pl.envelope[binY] = cfg.envelopeY;
  pl.envelope[binZ] = cfg.envelopeZ;
  return layerCreator.planeLayer(gctx, cfg.surfaces, cfg.binsY, cfg.binsZ,
                                 BinningValue::binX, pl, trafo);
}

std::pair<double, double> Acts::CuboidVolumeBuilder::binningRange(
    const GeometryContext& gctx,
    const Acts::CuboidVolumeBuilder::VolumeConfig& cfg) const {
  using namespace UnitLiterals;
  // Construct return value
  std::pair<double, double> minMax = std::make_pair(
      std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

  // Compute the min volume boundaries for computing the binning start
  // See
  // https://acts.readthedocs.io/en/latest/core/geometry.html#geometry-building
  // !! IMPORTANT !! The volume is assumed to be already rotated into the
  // telescope geometry
  Vector3 minVolumeBoundaries = cfg.position - 0.5 * cfg.length;
  Vector3 maxVolumeBoundaries = cfg.position + 0.5 * cfg.length;

  // Compute first the min-max from the layers

  for (const auto& layercfg : cfg.layerCfg) {
    // recreating the protolayer for each layer => slow, but only few sensors
    ProtoLayer pl{gctx, layercfg.surfaces};
    pl.envelope[binX] = layercfg.envelopeX;

    double surfacePosMin = pl.min(binX);
    double surfacePosMax = pl.max(binX);

    // Test if new extreme is found and set it
    if (surfacePosMin < minMax.first) {
      minMax.first = surfacePosMin;
    }
    if (surfacePosMax > minMax.second) {
      minMax.second = surfacePosMax;
    }
  }

  // Use the volume boundaries as limits for the binning
  minMax.first = std::min(minMax.first, minVolumeBoundaries(binX));
  minMax.second = std::max(minMax.second, maxVolumeBoundaries(binX));

  return minMax;
}

std::shared_ptr<Acts::TrackingVolume> Acts::CuboidVolumeBuilder::buildVolume(
    const GeometryContext& gctx,
    Acts::CuboidVolumeBuilder::VolumeConfig& cfg) const {
  // Build transformation
  Transform3 trafo(Transform3::Identity());
  trafo.translation() = cfg.position;
  // Set bounds
  auto bounds = std::make_shared<const CuboidVolumeBounds>(
      cfg.length.x() * 0.5, cfg.length.y() * 0.5, cfg.length.z() * 0.5);

  if (cfg.layerCfg.empty()) {
    // Build dummy layer if no layer is given (tmp solution)
    SurfaceConfig sCfg;
    sCfg.position = cfg.position;
    // Rotation of the surfaces: +pi/2 around axis y
    Vector3 xPos(0., 0., 1.);
    Vector3 yPos(0., 1., 0.);
    Vector3 zPos(-1., 0., 0.);
    sCfg.rotation.col(0) = xPos;
    sCfg.rotation.col(1) = yPos;
    sCfg.rotation.col(2) = zPos;
    // Bounds
    sCfg.rBounds = std::make_shared<const RectangleBounds>(
        RectangleBounds(cfg.length.y() * 0.5, cfg.length.z() * 0.5));

    LayerConfig lCfg;
    lCfg.surfaceCfg = {sCfg};

    cfg.layerCfg.push_back(lCfg);
  }

  // Gather the layers
  LayerVector layVec;
  if (cfg.layers.empty()) {
    cfg.layers.reserve(cfg.layerCfg.size());

    for (auto& layerCfg : cfg.layerCfg) {
      cfg.layers.push_back(buildLayer(gctx, layerCfg));
      layVec.push_back(cfg.layers.back());
    }
  } else {
    for (auto& lay : cfg.layers) {
      layVec.push_back(lay);
    }
  }

  // Build layer array
  std::pair<double, double> minMax = binningRange(gctx, cfg);
  LayerArrayCreator::Config lacCnf;
  LayerArrayCreator layArrCreator(
      lacCnf, getDefaultLogger("LayerArrayCreator", Logging::INFO));
  std::unique_ptr<const LayerArray> layArr(
      layArrCreator.layerArray(gctx, layVec, minMax.first, minMax.second,
                               BinningType::arbitrary, BinningValue::binX));

  // Build confined volumes
  if (cfg.trackingVolumes.empty()) {
    for (VolumeConfig vc : cfg.volumeCfg) {
      cfg.trackingVolumes.push_back(buildVolume(gctx, vc));
    }
  }

  std::shared_ptr<TrackingVolume> trackVolume;
  if (layVec.empty()) {
    // Build TrackingVolume
    trackVolume =
        TrackingVolume::create(trafo, bounds, cfg.volumeMaterial, nullptr,
                               nullptr, cfg.trackingVolumes, cfg.name);
  } else {
    // Build TrackingVolume
    trackVolume = TrackingVolume::create(trafo, bounds, cfg.volumeMaterial,
                                         std::move(layArr), nullptr,
                                         cfg.trackingVolumes, cfg.name);
  }
  return trackVolume;
}

Acts::MutableTrackingVolumePtr Acts::CuboidVolumeBuilder::trackingVolume(
    const GeometryContext& gctx, Acts::TrackingVolumePtr /*gctx*/,
    Acts::VolumeBoundsPtr /*bounds*/) const {
  // Build volumes
  std::vector<std::shared_ptr<TrackingVolume>> volumes;
  volumes.reserve(m_cfg.volumeCfg.size());
  for (VolumeConfig volCfg : m_cfg.volumeCfg) {
    volumes.push_back(buildVolume(gctx, volCfg));
  }

  // Sort the volumes vectors according to the center location, otherwise the
  // binning boundaries will fail
  std::sort(volumes.begin(), volumes.end(),
            [](const TrackingVolumePtr& lhs, const TrackingVolumePtr& rhs) {
              return lhs->center().x() < rhs->center().x();
            });

  // Glue volumes
  for (unsigned int i = 0; i < volumes.size() - 1; i++) {
    volumes[i + 1]->glueTrackingVolume(
        gctx, BoundarySurfaceFace::negativeFaceYZ, volumes[i].get(),
        BoundarySurfaceFace::positiveFaceYZ);
    volumes[i]->glueTrackingVolume(gctx, BoundarySurfaceFace::positiveFaceYZ,
                                   volumes[i + 1].get(),
                                   BoundarySurfaceFace::negativeFaceYZ);
  }

  // Translation
  Transform3 trafo(Transform3::Identity());
  trafo.translation() = m_cfg.position;

  // Size of the volume
  auto volume = std::make_shared<const CuboidVolumeBounds>(
      m_cfg.length.x() * 0.5, m_cfg.length.y() * 0.5, m_cfg.length.z() * 0.5);

  // Build vector of confined volumes
  std::vector<std::pair<TrackingVolumePtr, Vector3>> tapVec;
  tapVec.reserve(m_cfg.volumeCfg.size());
  for (auto& tVol : volumes) {
    tapVec.push_back(std::make_pair(tVol, tVol->center()));
  }

  // Set bin boundaries along binning
  std::vector<float> binBoundaries;
  binBoundaries.push_back(volumes[0]->center().x() -
                          m_cfg.volumeCfg[0].length.x() * 0.5);
  for (size_t i = 0; i < volumes.size(); i++) {
    binBoundaries.push_back(volumes[i]->center().x() +
                            m_cfg.volumeCfg[i].length.x() * 0.5);
  }

  // Build binning
  BinningData binData(BinningOption::open, BinningValue::binX, binBoundaries);
  auto bu = std::make_unique<const BinUtility>(binData);

  // Build TrackingVolume array
  std::shared_ptr<const TrackingVolumeArray> trVolArr(
      new BinnedArrayXD<TrackingVolumePtr>(tapVec, std::move(bu)));

  // Create world volume
  MutableTrackingVolumePtr mtvp(
      TrackingVolume::create(trafo, volume, trVolArr, "World"));

  return mtvp;
}
