// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <stdexcept>

namespace ActsDownstream {

struct GeometryModuleHandle {
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
};

inline std::unique_ptr<GeometryModuleHandle> buildTinyTrackingGeometryHandle() {
  using namespace Acts;
  using namespace Acts::UnitLiterals;

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const LayerArrayCreator>(
      lacConfig, getDefaultLogger("LayerArrayCreator", Logging::WARNING));

  TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator = std::make_shared<const TrackingVolumeArrayCreator>(
      tvacConfig,
      getDefaultLogger("TrackingVolumeArrayCreator", Logging::WARNING));

  CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper = std::make_shared<const CylinderVolumeHelper>(
      cvhConfig, getDefaultLogger("CylinderVolumeHelper", Logging::WARNING));

  PassiveLayerBuilder::Config layerBuilderConfig;
  layerBuilderConfig.layerIdentification = "DownstreamGeometryModule";
  layerBuilderConfig.centralLayerRadii = {10_mm};
  layerBuilderConfig.centralLayerHalflengthZ = {40_mm};
  layerBuilderConfig.centralLayerThickness = {1_mm};
  auto layerBuilder = std::make_shared<const PassiveLayerBuilder>(
      layerBuilderConfig,
      getDefaultLogger("DownstreamLayerBuilder", Logging::WARNING));

  CylinderVolumeBuilder::Config cvbConfig;
  cvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
  cvbConfig.volumeName = "DownstreamVolume";
  cvbConfig.layerBuilder = layerBuilder;
  cvbConfig.layerEnvelopeR = {1_mm, 1_mm};
  cvbConfig.buildToRadiusZero = true;
  auto centralVolumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
      cvbConfig, getDefaultLogger("DownstreamVolumeBuilder", Logging::WARNING));

  TrackingGeometryBuilder::Config tgbConfig;
  tgbConfig.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return centralVolumeBuilder->trackingVolume(context, inner);
      });
  tgbConfig.trackingVolumeHelper = cylinderVolumeHelper;

  TrackingGeometryBuilder tgBuilder(tgbConfig);
  auto trackingGeometry = tgBuilder.trackingGeometry(gctx);
  if (trackingGeometry == nullptr) {
    throw std::runtime_error(
        "Failed to build tracking geometry in downstream module");
  }

  auto handle = std::make_unique<GeometryModuleHandle>();
  handle->trackingGeometry = std::move(trackingGeometry);
  return handle;
}

inline void destroyTinyTrackingGeometryHandle(void* handle) {
  delete static_cast<GeometryModuleHandle*>(handle);
}

}  // namespace ActsDownstream
