// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/MockupSectorBuilder.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/PortalHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Plugins/Geant4/Geant4Converters.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

// a function for sorting the detector volumes from the innermost to the
// outermost

bool sortByRadialDistance(
    std::shared_ptr<Acts::Experimental::DetectorVolume>& det_vol1,
    std::shared_ptr<Acts::Experimental::DetectorVolume>& det_vol2) {
  return det_vol1->center().y() < det_vol2->center().y();
}

ActsExamples::MockupSectorBuilder::MockupSectorBuilder(
    const ActsExamples::MockupSectorBuilder::Config& config)
    : m_cfg() {
  setConfiguration(config);

  setWorld();
}

// initialize the m_cfg
void ActsExamples::MockupSectorBuilder::setConfiguration(
    const ActsExamples::MockupSectorBuilder::Config& config) {
  m_cfg = config;
}

// initiliaze the g4World
void ActsExamples::MockupSectorBuilder::setWorld() {
  ActsExamples::GdmlDetectorConstruction geo_gdml(m_cfg.gdml_path);
  g4World = geo_gdml.Construct();
}

std::shared_ptr<Acts::Experimental::DetectorVolume>
ActsExamples::MockupSectorBuilder::BuildChamber(
    const Acts::GeometryContext& gctx,
    const ActsExamples::MockupSectorBuilder::ChamberConfig& chamber_config) {
  if (g4World == nullptr) {
    throw std::invalid_argument(
        "MockupSector: No g4World initialized or number ");
    return nullptr;
  }

  // Geant4Detector Config creator with the g4world from the gdml file
  auto g4World_config = ActsExamples::Geant4::Geant4Detector::Config();
  g4World_config.name = "Chamber";
  g4World_config.g4World = g4World;

  // Get the sensitive and passive surfaces and pass to the g4World Config
  auto g4Sensitive =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          chamber_config.SensitiveNames);
  auto g4Passive =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          chamber_config.PassiveNames);

  auto g4SurfaceOptions = Acts::Geant4DetectorSurfaceFactory::Options();
  g4SurfaceOptions.sensitiveSurfaceSelector = g4Sensitive;
  g4SurfaceOptions.passiveSurfaceSelector = g4Passive;
  g4World_config.g4SurfaceOptions = g4SurfaceOptions;

  auto g4detector = ActsExamples::Geant4::Geant4Detector();

  auto [surfaces, detectorElements] =
      g4detector.convertGeant4Volumes(g4World_config);

  // The vector that holds the converted sensitive surfaces of the chamber
  std::vector<std::shared_ptr<Acts::Surface>> straw_surfaces = {};
  // The y and z coordinates of the straws' centers
  std::vector<float> positionY = {};
  std::vector<float> positionZ = {};
  std::vector<float> positionX = {};

  // Convert the physical volumes of the detector elements to straw surfaces
  for (auto& detectorElement : detectorElements) {
    auto context = Acts::GeometryContext();
    auto g4conv = Acts::Geant4PhysicalVolumeConverter();

    g4conv.forcedType = Acts::Surface::SurfaceType::Straw;
    auto g4conv_surf = g4conv.Geant4PhysicalVolumeConverter::surface(
        detectorElement->g4PhysicalVolume(),
        detectorElement->transform(context));

    straw_surfaces.push_back(g4conv_surf);
    positionY.push_back(g4conv_surf->center(context).y());
    positionZ.push_back(g4conv_surf->center(context).z());
    positionX.push_back(g4conv_surf->center(context).x());
  }

  // Create the bounds of the detector volumes
  float radius = straw_surfaces.front()->bounds().values()[0];
  auto minY = *min_element(positionY.begin(), positionY.end());
  auto maxY = *max_element(positionY.begin(), positionY.end());

  auto minZ = *min_element(positionZ.begin(), positionZ.end());
  auto maxZ = *max_element(positionZ.begin(), positionZ.end());

  auto minX = *min_element(positionX.begin(), positionX.end());
  auto maxX = *max_element(positionX.begin(), positionX.end());

  Acts::ActsScalar hx = straw_surfaces.front()->bounds().values()[1] + 10;
  Acts::ActsScalar hy = 0.5 * (std::abs(minY - maxY) + 2 * radius + 10.);
  Acts::ActsScalar hz = 0.5 * (std::abs(minZ - maxZ) + 2 * radius + 10.);

  auto detectorVolumeBounds =
      std::make_unique<Acts::CuboidVolumeBounds>(hx, hy, hz);

  Acts::Vector3 chamber_position = {(maxX + minX) / 2, (maxY + minY) / 2,
                                    (maxZ + minZ) / 2};

  // create the detector volume for the chamber
  auto detectorVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalGenerator(), gctx, chamber_config.name,
      Acts::Transform3(Acts::Translation3(chamber_position)),
      std::move(detectorVolumeBounds), straw_surfaces,
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>{},
      Acts::Experimental::allPortalsAndSurfaces());

  return detectorVolume;
}

std::shared_ptr<Acts::Experimental::DetectorVolume>
ActsExamples::MockupSectorBuilder::BuildSector(
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
        det_volumes,
    const Acts::GeometryContext& gctx) {
  if (m_cfg.NumberOfSectors > maxNumberOfSectors) {
    throw std::invalid_argument("MockupSector:Number of max sectors exceeded");
    return nullptr;
  }

  // sort the dtector volumes by their radial distance (from
  // innermost---->outermost)
  std::sort(det_volumes.begin(), det_volumes.end(), sortByRadialDistance);

  // calculate the inner and outer radius of the cylinderbounds
  auto rmin = det_volumes.front()->center().y() -
              2 * det_volumes.front()->volumeBounds().values()[1];
  auto rmax = det_volumes.back()->center().y() +
              3 * det_volumes.back()->volumeBounds().values()[1];
  auto halfZ = det_volumes.back()->volumeBounds().values()[2];

  // radius of the inner cylinder for the inner chambers
  auto rmin_inner = rmin;
  auto rmax_inner = det_volumes.front()->center().y() +
                    3 * det_volumes.front()->volumeBounds().values()[1];

  // radius of the middle cylinder for the middle chambers
  auto rmin_middle = det_volumes[1]->center().y() -
                     2 * det_volumes[1]->volumeBounds().values()[1];
  auto rmax_middle = det_volumes[1]->center().y() +
                     3 * det_volumes[1]->volumeBounds().values()[1];

  // radius of the outer cylinder
  auto rmin_outer = det_volumes.back()->center().y() -
                    2 * det_volumes.back()->volumeBounds().values()[1];
  auto rmax_outer = det_volumes.back()->center().y() +
                    3 * det_volumes.back()->volumeBounds().values()[1];

  // create two vectors in order to define the opening angle of each sector
  auto xA = det_volumes.back()->center().x() +
            det_volumes.back()->volumeBounds().values()[0];
  auto yA = det_volumes.back()->center().y() -
            det_volumes.back()->volumeBounds().values()[1];
  auto zA = det_volumes.back()->center().z();

  auto xB = det_volumes.back()->center().x() -
            det_volumes.back()->volumeBounds().values()[0];
  auto yB = det_volumes.back()->center().y() -
            det_volumes.back()->volumeBounds().values()[1];
  auto zB = det_volumes.back()->center().z();

  Acts::Vector3 pointA = {xA, yA, zA};
  Acts::Vector3 pointB = {xB, yB, zB};

  // calculate the phi angles of the vectors
  auto phiA = Acts::VectorHelpers::phi(pointA);
  auto phiB = Acts::VectorHelpers::phi(pointB);
  auto sector_angle = M_PI;

  auto hopening_angle = M_PI / m_cfg.NumberOfSectors;

  if (m_cfg.NumberOfSectors == 1) {
    hopening_angle = (phiB - phiA) / 2;
    sector_angle = hopening_angle;
  }

  auto detectorVolumeBounds = std::make_unique<Acts::CylinderVolumeBounds>(
      rmin, rmax, halfZ, sector_angle);
  auto inner_detectorVolumeBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(rmin_inner, rmax_inner,
                                                   halfZ, sector_angle);
  auto middle_detectorVolumeBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(rmin_middle, rmax_middle,
                                                   halfZ, sector_angle);
  auto outer_detectorVolumeBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(rmin_outer, rmax_outer,
                                                   halfZ, sector_angle);

  const Acts::Vector3 pos = {0., 0., 0.};

  // the transfom of the cylinder volume
  Acts::AngleAxis3 rotZ(M_PI / 2, Acts::Vector3(0., 0., 1));
  auto transform = Acts::Transform3(Acts::Translation3(pos));
  transform *= rotZ;

  // keep the inner, middle and outer chambers in order to constraint them in
  // cylinder bounds
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detectorVolumesOfSector = {};
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      inner_detectorVolumesOfSector = {};
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      middle_detectorVolumesOfSector = {};
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      outer_detectorVolumesOfSector = {};

  std::vector<std::shared_ptr<Acts::Surface>> shifted_surfaces = {};

  int itr = 0;

  for (int i = 0; i < m_cfg.NumberOfSectors; i++) {
    Acts::AngleAxis3 rotation(2 * i * hopening_angle,
                              Acts::Vector3(0., 0., 1.));

    itr = 0;

    for (auto& det_vol : det_volumes) {
      auto shift_vol =
          rotation * Acts::Transform3(Acts::Translation3(det_vol->center()));

      for (auto& det_surf : det_vol->surfaces()) {
        auto shift_surf = Acts::Transform3::Identity() * rotation;

        // create the shifted surfaces by creating copied surface objects
        auto radius = det_surf->bounds().values()[0];
        auto halflengthz = det_surf->bounds().values()[1];
        auto strawSurfaceObject = Acts::Surface::makeShared<Acts::StrawSurface>(
            det_surf->transform(Acts::GeometryContext()), radius, halflengthz);
        auto copiedTransformStrawSurface =
            Acts::Surface::makeShared<Acts::StrawSurface>(
                Acts::GeometryContext(), *strawSurfaceObject, shift_surf);
        shifted_surfaces.push_back(copiedTransformStrawSurface);
      }

      // create the bounds of the volumes of each chamber
      auto bounds = std::make_unique<Acts::CuboidVolumeBounds>(
          det_vol->volumeBounds().values()[0],
          det_vol->volumeBounds().values()[1],
          det_vol->volumeBounds().values()[2]);
      // create the shifted chamber
      auto detectorVolume_sec =
          Acts::Experimental::DetectorVolumeFactory::construct(
              Acts::Experimental::defaultPortalGenerator(), gctx,
              "test-cylinder", shift_vol, std::move(bounds), shifted_surfaces,
              std::vector<
                  std::shared_ptr<Acts::Experimental::DetectorVolume>>{},
              Acts::Experimental::allPortalsAndSurfaces());

      // detectorVolumesOfSector.push_back(detectorVolume_sec);

      if (itr == 0)
        inner_detectorVolumesOfSector.push_back(detectorVolume_sec);
      if (itr == 1)
        middle_detectorVolumesOfSector.push_back(detectorVolume_sec);
      if (itr == 2)
        outer_detectorVolumesOfSector.push_back(detectorVolume_sec);

      itr += 1;

      shifted_surfaces.clear();

    }  // end of detector volumes

  }  // end of number of sectors

  // creation of the cylinder volume of the inner chambers

  auto detectorVolume_inner =
      Acts::Experimental::DetectorVolumeFactory::construct(
          Acts::Experimental::defaultPortalGenerator(), gctx, "inner_cylinder",
          transform, std::move(inner_detectorVolumeBounds),
          std::vector<std::shared_ptr<Acts::Surface>>{},
          inner_detectorVolumesOfSector,
          Acts::Experimental::allPortalsAndSurfaces());

  auto detectorVolume_middle =
      Acts::Experimental::DetectorVolumeFactory::construct(
          Acts::Experimental::defaultPortalGenerator(), gctx, "middle_cylinder",
          transform, std::move(middle_detectorVolumeBounds),
          std::vector<std::shared_ptr<Acts::Surface>>{},
          middle_detectorVolumesOfSector,
          Acts::Experimental::allPortalsAndSurfaces());

  auto detectorVolume_outer =
      Acts::Experimental::DetectorVolumeFactory::construct(
          Acts::Experimental::defaultPortalGenerator(), gctx, "outer_cylinder",
          transform, std::move(outer_detectorVolumeBounds),
          std::vector<std::shared_ptr<Acts::Surface>>{},
          outer_detectorVolumesOfSector,
          Acts::Experimental::allPortalsAndSurfaces());

  detectorVolumesOfSector.push_back(detectorVolume_inner);
  detectorVolumesOfSector.push_back(detectorVolume_middle);
  detectorVolumesOfSector.push_back(detectorVolume_outer);

  std::cout << "volumes in sector" << detectorVolumesOfSector.size()
            << std::endl;

  // creation of the mother volume
  auto detectorVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalGenerator(), gctx, "sector", transform,
      std::move(detectorVolumeBounds),
      std::vector<std::shared_ptr<Acts::Surface>>{}, detectorVolumesOfSector,
      Acts::Experimental::allPortalsAndSurfaces());

  return detectorVolume;
}