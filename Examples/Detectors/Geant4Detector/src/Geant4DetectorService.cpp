// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/Geant4DetectorService.hpp"

#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/Detector.hpp"
#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/KDTreeTrackingGeometryBuilder.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include "G4Transform3D.hh"

ActsExamples::Geant4::Geant4DetectorService::Geant4DetectorService(
    const ActsExamples::Geant4::Geant4DetectorService::Config& cfg)
    : BareService(cfg.name, cfg.logLevel), m_cfg(cfg) {}

void ActsExamples::Geant4::Geant4DetectorService::startRun() {
  m_g4DetectorConstruction =
      std::make_unique<ActsExamples::GdmlDetectorConstruction>(m_cfg.gdmlFile);
  auto g4WorldVolume = m_g4DetectorConstruction->Construct();

  Acts::Geant4DetectorSurfaceFactory::Config g4DetElementConfig;
  Acts::Geant4DetectorSurfaceFactory g4DetElementFactory(g4DetElementConfig);

  Acts::Geant4DetectorSurfaceFactory::Cache g4DetElementCache;
  Acts::Geant4DetectorSurfaceFactory::Options g4DetElementOptions;

  G4Transform3D g4ToWorld;

  g4DetElementOptions.sensitiveSelector =
      Acts::Geant4PhysicalVolumeSelectors::generateNameSelector(
          m_cfg.sensitiveSelectionName);
  g4DetElementOptions.passiveSelector =
      Acts::Geant4PhysicalVolumeSelectors::generateNameSelector(
          m_cfg.passiveSelectionName);

  g4DetElementFactory.construct(g4DetElementCache, g4ToWorld, *g4WorldVolume,
                                g4DetElementOptions);

  ACTS_INFO("Found " << g4DetElementCache.matchedG4Volumes
                     << " matching  Geant4 Physical volumes.");

  ACTS_INFO("Found " << g4DetElementCache.convertedSurfaces
                     << " converted Geant4 Physical volumes.");

  ACTS_INFO("Found " << g4DetElementCache.convertedMaterials
                     << " converted Geant4 Material slabs.");

  if (m_cfg.buildTrackingGeometry and m_trackingGeometry == nullptr) {
    // Surface array creatorr
    auto surfaceArrayCreator =
        std::make_shared<const Acts::SurfaceArrayCreator>(
            Acts::SurfaceArrayCreator::Config(),
            Acts::getDefaultLogger("SurfaceArrayCreator", m_cfg.toolLogLevel));
    // Layer Creator
    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator = std::make_shared<Acts::LayerCreator>(
        lcConfig, Acts::getDefaultLogger("LayerCreator", m_cfg.toolLogLevel));
    // Layer array creator
    Acts::LayerArrayCreator::Config lacConfig;
    auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
        lacConfig,
        Acts::getDefaultLogger("LayerArrayCreator", m_cfg.toolLogLevel));
    // Tracking volume array creator
    Acts::TrackingVolumeArrayCreator::Config tvacConfig;
    auto tVolumeArrayCreator =
        std::make_shared<const Acts::TrackingVolumeArrayCreator>(
            tvacConfig, Acts::getDefaultLogger("TrackingVolumeArrayCreator",
                                               m_cfg.toolLogLevel));
    // configure the cylinder volume helper
    Acts::CylinderVolumeHelper::Config cvhConfig;
    cvhConfig.layerArrayCreator = layerArrayCreator;
    cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
    auto cylinderVolumeHelper =
        std::make_shared<const Acts::CylinderVolumeHelper>(
            cvhConfig,
            Acts::getDefaultLogger("CylinderVolumeHelper", m_cfg.toolLogLevel));

    // The KDT tracking geometry builder
    Acts::KDTreeTrackingGeometryBuilder::Config kdtgConfig;
    kdtgConfig.layerCreator = layerCreator;
    kdtgConfig.trackingVolumeHelper = cylinderVolumeHelper;
    // Reserve the right amount of surfaces
    std::vector<std::shared_ptr<Acts::Surface>> surfaces;
    kdtgConfig.surfaces.reserve(g4DetElementCache.sensitiveSurfaces.size() +
                                g4DetElementCache.passiveSurfaces.size());
    // Add the sensitive surfaces
    for (const auto& e : g4DetElementCache.sensitiveSurfaces) {
      kdtgConfig.surfaces.push_back(
          std::get<std::shared_ptr<Acts::Surface>>(e));
    }
    // Add the passive surfaces
    kdtgConfig.surfaces.insert(kdtgConfig.surfaces.end(),
                               g4DetElementCache.passiveSurfaces.begin(),
                               g4DetElementCache.passiveSurfaces.end());

    // Assign the proto detector
    kdtgConfig.protoDetector = m_cfg.protoDetector;

    // Make the builder
    auto kdtTrackingGeometryBuilder = Acts::KDTreeTrackingGeometryBuilder(
        kdtgConfig, Acts::getDefaultLogger("KDTreeTrackingGeometryBuilder",
                                           m_cfg.toolLogLevel));

    Acts::GeometryContext tContext;
    m_trackingGeometry = kdtTrackingGeometryBuilder.trackingGeometry(tContext);
  }
}
