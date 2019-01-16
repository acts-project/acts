// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Tools/CuboidVolumeBuilder.hpp"

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Detector/detail/DefaultDetectorElementBase.hpp"
#include "Acts/Layers/Layer.hpp"
#include "Acts/Layers/PlaneLayer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Tools/LayerCreator.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/CuboidVolumeBounds.hpp"

std::shared_ptr<const Acts::PlaneSurface>
Acts::CuboidVolumeBuilder::buildSurface(
    const Acts::CuboidVolumeBuilder::SurfaceConfig& cfg) const
{
  std::shared_ptr<PlaneSurface> surface;

  // Build transformation
  Transform3D trafo(Transform3D::Identity() * cfg.rotation);
  trafo.translation() = cfg.position;

  // Create and store surface
  if (cfg.detElementConstructor) {
    surface = Surface::makeShared<PlaneSurface>(
        cfg.rBounds,
        *(cfg.detElementConstructor(std::make_shared<const Transform3D>(trafo),
                                    cfg.rBounds,
                                    cfg.thickness)));
  } else {
    surface = Surface::makeShared<PlaneSurface>(
        std::make_shared<const Transform3D>(trafo), cfg.rBounds);
  }
  surface->setAssociatedMaterial(cfg.surMat);
  return surface;
}

std::shared_ptr<const Acts::Layer>
Acts::CuboidVolumeBuilder::buildLayer(
    Acts::CuboidVolumeBuilder::LayerConfig& cfg) const
{
  // Build the surface
  if (cfg.surface == nullptr) {
    cfg.surface = buildSurface(cfg.surfaceCfg);
  }
  // Build transformation centered at the surface position
  Transform3D trafo(Transform3D::Identity() * cfg.surfaceCfg.rotation);
  trafo.translation() = cfg.surfaceCfg.position;

  LayerCreator::Config lCfg;
  lCfg.surfaceArrayCreator = std::make_shared<const SurfaceArrayCreator>();
  LayerCreator layerCreator(lCfg);

  return layerCreator.planeLayer({cfg.surface},
                                 cfg.binsY,
                                 cfg.binsZ,
                                 BinningValue::binX,
                                 boost::none,
                                 std::make_shared<const Transform3D>(trafo));
}

std::pair<double, double>
Acts::CuboidVolumeBuilder::binningRange(
    const Acts::CuboidVolumeBuilder::VolumeConfig& cfg) const
{
  // Construct return value
  std::pair<double, double> minMax = std::make_pair(
      std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
  for (const auto& layercfg : cfg.layerCfg) {
    // Test if new extreme is found and set it
    if (layercfg.surfaceCfg.position.x() - 1. * units::_um < minMax.first) {
      minMax.first = layercfg.surfaceCfg.position.x() - 1. * units::_um;
    }
    if (layercfg.surfaceCfg.position.x() + 1. * units::_um > minMax.second) {
      minMax.second = layercfg.surfaceCfg.position.x() + 1. * units::_um;
    }
  }
  return minMax;
}

std::shared_ptr<Acts::TrackingVolume>
Acts::CuboidVolumeBuilder::buildVolume(
    Acts::CuboidVolumeBuilder::VolumeConfig& cfg) const
{
  // Build transformation
  Transform3D trafo(Transform3D::Identity());
  trafo.translation() = cfg.position;
  // Set bounds
  auto bounds = std::make_shared<const CuboidVolumeBounds>(
      cfg.length.x() * 0.5, cfg.length.y() * 0.5, cfg.length.z() * 0.5);

  if (cfg.layerCfg.empty()) {
    // Build dummy layer if no layer is given (tmp solution)
    SurfaceConfig sCfg;
    sCfg.position = cfg.position;
    // Rotation of the surfaces: +pi/2 around axis y
    Vector3D xPos(0., 0., 1.);
    Vector3D yPos(0., 1., 0.);
    Vector3D zPos(-1., 0., 0.);
    sCfg.rotation.col(0) = xPos;
    sCfg.rotation.col(1) = yPos;
    sCfg.rotation.col(2) = zPos;
    // Bounds
    sCfg.rBounds = std::make_shared<const RectangleBounds>(
        RectangleBounds(cfg.length.y() * 0.5, cfg.length.z() * 0.5));

    LayerConfig lCfg;
    lCfg.surfaceCfg = sCfg;

    cfg.layerCfg.push_back(lCfg);
  }

  // Gather the layers
  LayerVector layVec;
  if (cfg.layers.empty()) {
    cfg.layers.reserve(cfg.layerCfg.size());

    for (auto& layerCfg : cfg.layerCfg) {
      cfg.layers.push_back(buildLayer(layerCfg));
      layVec.push_back(cfg.layers.back());
    }
  } else {
    for (auto& lay : cfg.layers) {
      layVec.push_back(lay);
    }
  }

  // Build layer array
  std::pair<double, double> minMax = binningRange(cfg);
  LayerArrayCreator layArrCreator(
      getDefaultLogger("LayerArrayCreator", Logging::INFO));
  std::unique_ptr<const LayerArray> layArr(
      layArrCreator.layerArray(layVec,
                               minMax.first,
                               minMax.second,
                               BinningType::arbitrary,
                               BinningValue::binX));

  // Build TrackingVolume
  auto trackVolume
      = TrackingVolume::create(std::make_shared<const Transform3D>(trafo),
                               bounds,
                               cfg.material,
                               std::move(layArr),
                               layVec,
                               {},
                               {},
                               cfg.name);
  trackVolume->sign(GeometrySignature::Global);

  return trackVolume;
}

Acts::MutableTrackingVolumePtr Acts::CuboidVolumeBuilder::trackingVolume(
    Acts::TrackingVolumePtr /*unused*/,
    Acts::VolumeBoundsPtr /*unused*/) const
{
  // Build volumes
  std::vector<std::shared_ptr<TrackingVolume>> volumes;
  volumes.reserve(m_cfg.volumeCfg.size());
  for (VolumeConfig volCfg : m_cfg.volumeCfg) {
    volumes.push_back(buildVolume(volCfg));
  }

  // Glue volumes
  for (unsigned int i = 0; i < volumes.size() - 1; i++) {
    volumes[i + 1]->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
                                       volumes[i],
                                       BoundarySurfaceFace::positiveFaceYZ);
    volumes[i]->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
                                   volumes[i + 1],
                                   BoundarySurfaceFace::negativeFaceYZ);
  }

  // Translation
  Transform3D trafo(Transform3D::Identity());
  trafo.translation() = m_cfg.position;

  // Size of the volume
  auto volume = std::make_shared<const CuboidVolumeBounds>(
      m_cfg.length.x() * 0.5, m_cfg.length.y() * 0.5, m_cfg.length.z() * 0.5);

  // Build vector of confined volumes
  std::vector<std::pair<TrackingVolumePtr, Vector3D>> tapVec;
  tapVec.reserve(m_cfg.volumeCfg.size());
  for (auto& tVol : volumes) {
    tapVec.push_back(std::make_pair(tVol, tVol->center()));
  }

  // Set bin boundaries along binning
  std::vector<float> binBoundaries;
  binBoundaries.push_back(volumes[0]->center().x()
                          - m_cfg.volumeCfg[0].length.x() * 0.5);
  for (size_t i = 0; i < volumes.size(); i++) {
    binBoundaries.push_back(volumes[i]->center().x()
                            + m_cfg.volumeCfg[i].length.x() * 0.5);
  }

  // Build binning
  BinningData binData(BinningOption::open, BinningValue::binX, binBoundaries);
  auto        bu = std::make_unique<const BinUtility>(binData);

  // Build TrackingVolume array
  std::shared_ptr<const TrackingVolumeArray> trVolArr(
      new BinnedArrayXD<TrackingVolumePtr>(tapVec, std::move(bu)));

  // Create world volume
  MutableTrackingVolumePtr mtvp(TrackingVolume::create(
      std::make_shared<const Transform3D>(trafo), volume, trVolArr, "World"));

  mtvp->sign(GeometrySignature::Global);
  return mtvp;
}
