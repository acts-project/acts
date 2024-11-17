// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"

#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/KDTreeTrackingGeometryBuilder.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"

#include <ostream>
#include <stdexcept>

#include "G4Transform3D.hh"
#include "G4VPhysicalVolume.hh"

namespace ActsExamples::Geant4 {

Geant4Detector::Geant4Detector(const Geant4Detector::Config& cfg)
    : DetectorCommons::Detector(
          Acts::getDefaultLogger("Geant4Detector", cfg.logLevel)),
      m_cfg(cfg) {}

const Geant4Detector::DetectorElements& Geant4Detector::detectorElements()
    const {
  return m_detectorElements;
}

void Geant4Detector::drop() {
  Detector::drop();

  m_detectorElements.clear();
}

void Geant4Detector::buildTrackingGeometry() {
  std::tie(m_trackingGeometry, m_contextDecorators, m_detectorElements) =
      constructTrackingGeometry();
}

void Geant4Detector::buildDetector() {
  std::tie(m_detector, m_contextDecorators, m_detectorElements) =
      constructDetector();
}

std::tuple<Geant4Detector::DetectorPtr, Geant4Detector::ContextDecorators,
           Geant4Detector::DetectorElements>
Geant4Detector::constructDetector() const {
  if (m_cfg.g4World == nullptr) {
    throw std::invalid_argument(
        "Geant4Detector: no world Geant4 volume provided");
  }

  ACTS_INFO("Building an Acts::Detector called '"
            << m_cfg.name << "' from the Geant4PhysVolume '"
            << m_cfg.g4World->GetName() << "'");

  DetectorPtr detector = nullptr;
  ContextDecorators decorators = {};

  auto [surfaces, elements] = convertGeant4Volumes();

  return {std::move(detector), std::move(decorators), std::move(elements)};
}

std::tuple<Geant4Detector::TrackingGeometryPtr,
           Geant4Detector::ContextDecorators, Geant4Detector::DetectorElements>
Geant4Detector::constructTrackingGeometry() const {
  if (m_cfg.g4World == nullptr) {
    throw std::invalid_argument(
        "Geant4Detector: no world Geant4 volume provided");
  }

  ACTS_INFO("Building an Acts::TrackingGeometry called '"
            << m_cfg.name << "' from the Geant4PhysVolume '"
            << m_cfg.g4World->GetName() << "'");

  ContextDecorators decorators = {};

  auto [surfaces, elements] = convertGeant4Volumes();

  // Surface array creator
  auto surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(
      Acts::SurfaceArrayCreator::Config(),
      logger().clone("SurfaceArrayCreator"));
  // Layer Creator
  Acts::LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<Acts::LayerCreator>(
      lcConfig, logger().clone("LayerCreator"));
  // Layer array creator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig, logger().clone("LayerArrayCreator"));
  // Tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig, logger().clone("TrackingVolumeArrayCreator"));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig, logger().clone("CylinderVolumeHelper"));

  // Configure the tracking geometry builder, copy the surfaces in
  Acts::KDTreeTrackingGeometryBuilder::Config kdtCfg;
  kdtCfg.surfaces = surfaces;
  kdtCfg.layerCreator = layerCreator;
  kdtCfg.trackingVolumeHelper = cylinderVolumeHelper;
  kdtCfg.protoDetector = m_cfg.protoDetector;
  kdtCfg.geometryIdentifierHook = m_cfg.geometryIdentifierHook;

  // The KDT tracking geometry builder
  auto kdtBuilder = Acts::KDTreeTrackingGeometryBuilder(
      kdtCfg, logger().clone("KDTreeTrackingGeometryBuilder"));

  Acts::GeometryContext tContext;
  TrackingGeometryPtr trackingGeometry = kdtBuilder.trackingGeometry(tContext);

  return {std::move(trackingGeometry), std::move(decorators),
          std::move(elements)};
}

std::tuple<Geant4Detector::Surfaces, Geant4Detector::DetectorElements>
Geant4Detector::convertGeant4Volumes() const {
  // Generate the surface cache
  Acts::Geant4DetectorSurfaceFactory::Cache g4SurfaceCache;
  G4Transform3D g4ToWorld;

  Acts::Geant4DetectorSurfaceFactory{}.construct(
      g4SurfaceCache, g4ToWorld, *m_cfg.g4World, m_cfg.g4SurfaceOptions);

  ACTS_INFO("Found " << g4SurfaceCache.matchedG4Volumes
                     << " matching  Geant4 Physical volumes.");
  ACTS_INFO("Found " << g4SurfaceCache.sensitiveSurfaces.size()
                     << " converted sensitive Geant4 Physical volumes.");
  ACTS_INFO("Found " << g4SurfaceCache.passiveSurfaces.size()
                     << " converted passive Geant4 Physical volumes.");
  ACTS_INFO("Found " << g4SurfaceCache.convertedMaterials
                     << " converted Geant4 Material slabs.");

  Surfaces surfaces = {};
  DetectorElements elements = {};

  // Reserve the right amount of surfaces
  surfaces.reserve(g4SurfaceCache.sensitiveSurfaces.size() +
                   g4SurfaceCache.passiveSurfaces.size());
  elements.reserve(g4SurfaceCache.sensitiveSurfaces.size());

  // Add the sensitive surfaces
  for (const auto& [e, s] : g4SurfaceCache.sensitiveSurfaces) {
    elements.push_back(e);
    surfaces.push_back(s);
  }
  // Add the passive surfaces
  surfaces.insert(surfaces.end(), g4SurfaceCache.passiveSurfaces.begin(),
                  g4SurfaceCache.passiveSurfaces.end());

  return {std::move(surfaces), std::move(elements)};
}

}  // namespace ActsExamples::Geant4
