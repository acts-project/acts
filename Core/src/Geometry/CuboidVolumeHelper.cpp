// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CuboidVolumeHelper.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/ILayerArrayCreator.hpp"
#include "Acts/Geometry/ITrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <tuple>

namespace {
/// @brief Helper method to create bounds and transform
///
/// @param dimension the dimensions parameters
///
/// @return a tuple of a bounds object and a transform
std::tuple<std::shared_ptr<Acts::CuboidVolumeBounds>, Acts::Transform3>
boundsAndTransform(const Acts::Extent& dimension) {
  // Create bounds
  std::array<Acts::BinningValue, 3u> carts = {Acts::binX, Acts::binY,
                                              Acts::binZ};
  std::array<double, 3> boundValues;
  Acts::Vector3 position(0., 0., 0.);
  for (const auto& bv : carts) {
    boundValues[bv] = dimension.span(bv) * 0.5;
    position[bv] = dimension.medium(bv);
  }
  auto volumeBounds = std::make_shared<Acts::CuboidVolumeBounds>(boundValues);
  Acts::Translation3 translation(position);
  auto transform = Acts::Transform3(translation);
  return std::tie(volumeBounds, transform);
}

/// @brief Create a bunch of layers
///
/// @param positions the positions where to create them
/// @param rotation the rotation where to create them
/// @param layerBounds the rectangular bounds
/// @param thickness the thickness of the layer
///
/// @return a bunch of layers
std::vector<std::shared_ptr<const Acts::Layer>> createLayers(
    const std::vector<Acts::Vector3> positions,
    const Acts::RotationMatrix3 rotation,
    std::shared_ptr<const Acts::RectangleBounds> layerBounds,
    Acts::ActsScalar thickness) {
  std::vector<std::shared_ptr<const Acts::Layer>> layers = {};
  layers.reserve(positions.size());
  for (const auto& p : positions) {
    Acts::Transform3 transform = Acts::Transform3::Identity();
    transform.pretranslate(p);
    transform.prerotate(rotation);
    layers.push_back(
        Acts::PlaneLayer::create(transform, layerBounds, nullptr, thickness));
  }
  return layers;
}

/// @brief Helper for layer rotation and reference position
///
/// @param dimension the volume extent
/// @param envelope the enveloope
/// @param bValue the binning value
///
/// @return a tuple of a rotation matrix and a reference position
std::tuple<Acts::RotationMatrix3, std::shared_ptr<Acts::RectangleBounds>>
layerRotationAndBounds(const Acts::Extent& dimension, Acts::ActsScalar envelope,
                       Acts::BinningValue bValue) {
  std::shared_ptr<Acts::RectangleBounds> bounds = nullptr;
  auto rotation = Acts::Transform3::Identity().rotation();
  if (bValue == Acts::binX) {
    bounds = std::make_shared<Acts::RectangleBounds>(
        dimension.span(Acts::binY) * 0.5 - envelope,
        dimension.span(Acts::binZ) * 0.5 - envelope);
    Acts::Vector3 colX = rotation.col(0);
    Acts::Vector3 colY = rotation.col(1);
    Acts::Vector3 colZ = rotation.col(2);
    rotation.col(0) = colY;
    rotation.col(1) = colZ;
    rotation.col(2) = colX;
  } else if (bValue == Acts::binY) {
    bounds = std::make_shared<Acts::RectangleBounds>(
        dimension.span(Acts::binZ) * 0.5 - envelope,
        dimension.span(Acts::binX) * 0.5 - envelope);
    Acts::Vector3 colX = rotation.col(0);
    Acts::Vector3 colY = rotation.col(1);
    Acts::Vector3 colZ = rotation.col(2);
    rotation.col(0) = colZ;
    rotation.col(1) = colX;
    rotation.col(2) = colY;
  } else {
    bounds = std::make_shared<Acts::RectangleBounds>(
        dimension.span(Acts::binX) * 0.5 - envelope,
        dimension.span(Acts::binY) * 0.5 - envelope);
  }
  return std::tie(rotation, bounds);
}

}  // namespace

Acts::CuboidVolumeHelper::CuboidVolumeHelper(
    const Acts::CuboidVolumeHelper::Config& cvhConfig,
    std::unique_ptr<const Logger> logger)
    : Acts::ITrackingVolumeHelper(), m_cfg(), m_logger(std::move(logger)) {
  setConfiguration(cvhConfig);
  // Throw exception if misconfigured
  if (static_cast<int>(m_cfg.bValue) > 2) {
    throw std::invalid_argument(
        "CuboidVolumeHelper: binning must be x, y or z");
  }
}

// Configuration
void Acts::CuboidVolumeHelper::setConfiguration(
    const Acts::CuboidVolumeHelper::Config& cvhConfig) {
  // copy the configuration
  m_cfg = cvhConfig;
  // Throw exception if misconfigured
  if (static_cast<int>(m_cfg.bValue) > 2) {
    throw std::invalid_argument(
        "CuboidVolumeHelper: binning must be x, y or z");
  }
}

void Acts::CuboidVolumeHelper::setLogger(
    std::unique_ptr<const Logger> newLogger) {
  m_logger = std::move(newLogger);
}

std::shared_ptr<Acts::TrackingVolume>
Acts::CuboidVolumeHelper::createTrackingVolume(
    const GeometryContext& gctx, const LayerVector& layers,
    std::shared_ptr<const IVolumeMaterial> volumeMaterial,
    std::shared_ptr<const VolumeBounds> volumeBounds,
    MutableTrackingVolumeVector mtvVector, const Transform3& transform,
    const std::string& volumeName, BinningType bType) const {
  std::unique_ptr<const LayerArray> layerArray = nullptr;
  if (not layers.empty()) {
    // Set min / max
    ActsScalar bPos =
        VectorHelpers::cast(transform.translation(), m_cfg.bValue);
    auto boundValues = volumeBounds->values();
    ActsScalar min = bPos - boundValues[m_cfg.bValue];
    ActsScalar max = bPos + boundValues[m_cfg.bValue];

    layerArray = m_cfg.layerArrayCreator->layerArray(gctx, layers, min, max,
                                                     bType, m_cfg.bValue);
  }

  // finally create the TrackingVolume
  auto tVolume = TrackingVolume::create(transform, volumeBounds, volumeMaterial,
                                        std::move(layerArray), nullptr,
                                        mtvVector, volumeName);
  // Screen output
  ACTS_VERBOSE(
      "Created cubic volume at position :" << toString(tVolume->center()));
  ACTS_VERBOSE("   created bounds : " << tVolume->volumeBounds());
  // Return the constructed TrackingVolume
  return tVolume;
}

std::shared_ptr<Acts::TrackingVolume>
Acts::CuboidVolumeHelper::createTrackingVolume(
    const GeometryContext& gctx, const LayerVector& layers,
    MutableTrackingVolumeVector mtvVector,
    std::shared_ptr<const IVolumeMaterial> volumeMaterial,
    const Extent& dimension, const std::string& volumeName,
    BinningType bType) const {
  // Check & create bounds and transform
  if (not dimension.constrains(binX) and not dimension.constrains(binY) and
      not dimension.constrains(binZ)) {
    throw std::invalid_argument(
        "CuboidVolumeHelper: x, y and z need to be constrained.");
  }
  auto [volumeBounds, transform] = boundsAndTransform(dimension);
  // Return the volume
  return createTrackingVolume(gctx, layers, volumeMaterial, volumeBounds,
                              mtvVector, transform, volumeName, bType);
}

std::shared_ptr<Acts::TrackingVolume>
Acts::CuboidVolumeHelper::createGapTrackingVolume(
    const GeometryContext& /*gctx*/, MutableTrackingVolumeVector& mtvVector,
    std::shared_ptr<const IVolumeMaterial> volumeMaterial,
    const Extent& dimension, unsigned int materialLayers,
    Surface::SurfaceType /*layerType*/, const std::string& volumeName) const {
  // Check & create bounds and transform
  if (not dimension.constrains(binX) and not dimension.constrains(binY) and
      not dimension.constrains(binZ)) {
    throw std::invalid_argument(
        "CuboidVolumeHelper: x, y and z need to be constrained.");
  }
  auto [volumeBounds, transform] = boundsAndTransform(dimension);

  std::unique_ptr<const LayerArray> layerArray = nullptr;
  if (materialLayers > 0) {
    auto [layerRotation, layerBounds] =
        layerRotationAndBounds(dimension, m_cfg.envelope, m_cfg.bValue);
    std::vector<Acts::Vector3> lPs;
    lPs.reserve(materialLayers);
    Acts::ActsScalar lStep =
        dimension.span(m_cfg.bValue) / (materialLayers + 1);
    for (unsigned int ilp = 1; ilp <= materialLayers; ++ilp) {
      Acts::Vector3 clp = transform.translation();
      clp[m_cfg.bValue] = dimension.min(m_cfg.bValue) + ilp * lStep;
      lPs.push_back(clp);
    }
    auto layers = createLayers(lPs, layerRotation, layerBounds,
                               m_cfg.passiveLayerThickness);
  }
  // Finally create the TrackingVolume
  return TrackingVolume::create(transform, volumeBounds, volumeMaterial,
                                std::move(layerArray), nullptr, mtvVector,
                                volumeName);
}

std::shared_ptr<Acts::TrackingVolume>
Acts::CuboidVolumeHelper::createGapTrackingVolume(
    const GeometryContext& gctx, MutableTrackingVolumeVector& mtvVector,
    std::shared_ptr<const IVolumeMaterial> volumeMaterial,
    const Extent& dimension, const std::vector<double>& layerPositions,
    Surface::SurfaceType /*layerType*/, const std::string& volumeName,
    BinningType bType) const {
  // Check & create bounds and transform
  if (not dimension.constrains(binX) and not dimension.constrains(binY) and
      not dimension.constrains(binZ)) {
    throw std::invalid_argument(
        "CuboidVolumeHelper: x, y and z need to be constrained.");
  }
  auto [volumeBounds, transform] = boundsAndTransform(dimension);

  std::unique_ptr<const LayerArray> layerArray = nullptr;
  if (not layerPositions.empty()) {
    auto [layerRotation, layerBounds] =
        layerRotationAndBounds(dimension, m_cfg.envelope, m_cfg.bValue);
    std::vector<Acts::Vector3> lPs;
    lPs.reserve(layerPositions.size());
    for (const auto lp : layerPositions) {
      Acts::Vector3 clp = transform.translation();
      clp[m_cfg.bValue] = lp;
      lPs.push_back(clp);
    }
    auto layers = createLayers(lPs, layerRotation, layerBounds,
                               m_cfg.passiveLayerThickness);

    // Set min / max
    ActsScalar bPos =
        VectorHelpers::cast(transform.translation(), m_cfg.bValue);
    auto boundValues = volumeBounds->values();
    ActsScalar min = bPos - boundValues[m_cfg.bValue];
    ActsScalar max = bPos + boundValues[m_cfg.bValue];
    layerArray = m_cfg.layerArrayCreator->layerArray(gctx, layers, min, max,
                                                     bType, m_cfg.bValue);
  }

  // Finally create the TrackingVolume
  return TrackingVolume::create(transform, volumeBounds, volumeMaterial,
                                std::move(layerArray), nullptr, mtvVector,
                                volumeName);
}

std::shared_ptr<Acts::TrackingVolume>
Acts::CuboidVolumeHelper::createContainerTrackingVolume(
    const GeometryContext& /*gctx*/,
    const TrackingVolumeVector& volumes) const {
  // Find out the binning direction through last and first volume
  const Vector3 firstPosition = volumes.front()->center();
  const Vector3 lastPosition = volumes.back()->center();
  const Vector3 diffPositions = lastPosition - firstPosition;

  // Guarantee the binnning
  BinningValue bValue = binZ;
  std::array<BinningValue, 2u> bValuesT = {binX, binY};
  std::array<BoundarySurfaceFace, 2u> glueFaces = {negativeFaceXY,
                                                   positiveFaceXY};

  if (diffPositions.x() > diffPositions.y() and
      diffPositions.x() > diffPositions.z()) {
    bValue = binX;
    bValuesT = {binY, binZ};
    glueFaces = {negativeFaceYZ, positiveFaceYZ};
  } else if (diffPositions.y() > diffPositions.x() and
             diffPositions.y() > diffPositions.z()) {
    bValue = binY;
    bValuesT = {binZ, binX};
    glueFaces = {negativeFaceZX, positiveFaceZX};
  }

  // For the binning
  std::vector<float> binSteps = {};
  ActsScalar binStart = firstPosition[bValue];

  std::array<ActsScalar, 3u> maxHalfs = {
      0.,
      0.,
      0.,
  };

  std::shared_ptr<TrackingVolume> lastVolume = nullptr;
  ActsScalar lastArea = 0.;

  std::vector<std::pair<std::shared_ptr<const TrackingVolume>, Vector3>> taps;
  for (auto [iv, v] : enumerate(volumes)) {
    // Get the overall bound values
    auto boundValues = v->volumeBounds().values();
    if (v->volumeBounds().type() != VolumeBounds::BoundsType::eCuboid) {
      throw std::runtime_error(
          "CuboidVolumeHelper: volume bounds must be cubic.");
    }

    // First volume, set the start value
    if (iv == 0) {
      binStart -= boundValues[bValue];
      binSteps.push_back(static_cast<float>(binStart));
    }

    for (auto [ibv, bv] : enumerate(boundValues)) {
      maxHalfs[ibv] = std::max(maxHalfs[ibv], bv);
    }
    // Glue them
    auto thisVolume = std::const_pointer_cast<TrackingVolume>(v);
    ActsScalar thisArea =
        4 * (boundValues[bValuesT[0u]] * boundValues[bValuesT[1u]]);
    if (lastVolume != nullptr) {
      if (thisArea > lastArea) {
        auto keeper = thisVolume->boundarySurfaces()[glueFaces[0u]];
        lastVolume->updateBoundarySurface(glueFaces[1u], keeper);
      } else {
        auto keeper = lastVolume->boundarySurfaces()[glueFaces[1u]];
        thisVolume->updateBoundarySurface(glueFaces[0u], keeper);
      }
    }
    // Remember the volumes
    lastVolume = thisVolume;
    lastArea = thisArea;
    binSteps.push_back(static_cast<float>(binStart + 2 * boundValues[bValue]));
    binStart += 2 * boundValues[bValue];
    // For the filling
    taps.push_back({v, v->center()});
  }

  // Evaluate the position
  Vector3 cPosition(0., 0., 0.);
  cPosition[bValue] = 0.5 * (binSteps.front() + binSteps.back());
  Translation3 cTranslation(cPosition);
  // Set the max half of the binning direction
  maxHalfs[bValue] = 0.5 * (binSteps.back() - binSteps.front());

  // Create new bounds
  auto cBounds = std::make_shared<CuboidVolumeBounds>(maxHalfs);
  auto cBinning = std::make_unique<BinUtility>(binSteps, open, bValue);
  using ConfinedVolumes = BinnedArrayXD<std::shared_ptr<const TrackingVolume>>;
  auto cVolumes = std::make_shared<ConfinedVolumes>(taps, std::move(cBinning));

  return TrackingVolume::create(Transform3(cTranslation), cBounds, cVolumes,
                                "Container");
}
