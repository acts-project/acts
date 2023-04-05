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
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <utility>  


ActsExamples::MockupSectorBuilder::MockupSectorBuilder(
    const ActsExamples::MockupSectorBuilder::Config& config){

        m_cfg = config;
       ActsExamples::GdmlDetectorConstruction geo_gdml(m_cfg.gdml_path);
       g4World = geo_gdml.Construct();
}


std::shared_ptr<Acts::Experimental::DetectorVolume>
ActsExamples::MockupSectorBuilder::buildChamber(
    const Acts::GeometryContext& gctx,
    const ActsExamples::MockupSectorBuilder::ChamberConfig& chamber_config) {
  if (g4World == nullptr) {
    throw std::invalid_argument(
        "MockupSector: No g4World initialized");
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

  auto [detector, surfaces, detectorElements] =
      g4detector.constructDetector(g4World_config, Acts::getDummyLogger());

  // The vector that holds the converted sensitive surfaces of the chamber
  std::vector<std::shared_ptr<Acts::Surface>> straw_surfaces = {};

   //auto min_max = std::pair<float, float>(std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
  auto max_limit = std::numeric_limits<float>::max();
  auto min_limit = -std::numeric_limits<float>::max();

   std::array<std::pair <float,float>,3> min_max;
   std::fill(min_max.begin(), min_max.end(), 
            std::make_pair<float, float>(std::move(max_limit), std::move(min_limit)));


  // Convert the physical volumes of the detector elements to straw surfaces
  for (auto& detectorElement : detectorElements) {
    auto context = Acts::GeometryContext();
    auto g4conv = Acts::Geant4PhysicalVolumeConverter();

    g4conv.forcedType = Acts::Surface::SurfaceType::Straw;
    auto g4conv_surf = g4conv.Geant4PhysicalVolumeConverter::surface(
        detectorElement->g4PhysicalVolume(),
        detectorElement->transform(context));

    straw_surfaces.push_back(g4conv_surf);

    min_max[0].first = std::min(min_max[0].first, (float) g4conv_surf->center(context).x());
    min_max[0].second = std::max(min_max[0].second, (float) g4conv_surf->center(context).x());

    min_max[1].first = std::min(min_max[1].first, (float) g4conv_surf->center(context).y());
    min_max[1].second = std::max(min_max[1].second, (float) g4conv_surf->center(context).y());

    min_max[2].first = std::min(min_max[2].first, (float) g4conv_surf->center(context).z());
    min_max[2].second = std::max(min_max[2].second, (float) g4conv_surf->center(context).z());

  }

  // Create the bounds of the detector volumes
  float radius = straw_surfaces.front()->bounds().values()[0];
  
  Acts::Vector3 minValues = {min_max[0].first, min_max[1].first, min_max[2].first};
  Acts::Vector3 maxValues = {min_max[0].second, min_max[1].second, min_max[2].second};


  Acts::ActsScalar hx = straw_surfaces.front()->bounds().values()[1] + m_cfg.toleranceOverlap;
  Acts::ActsScalar hy = 0.5 * (std::abs(minValues.y() - maxValues.y()) + 2 * radius + m_cfg.toleranceOverlap);
  Acts::ActsScalar hz = 0.5 * (std::abs(minValues.z() - maxValues.z()) + 2 * radius + m_cfg.toleranceOverlap);

  auto detectorVolumeBounds =
      std::make_unique<Acts::CuboidVolumeBounds>(hx, hy, hz);

  Acts::Vector3 chamber_position = {(maxValues.x() + minValues.x()) / 2, (maxValues.y() + minValues.y()) / 2,
                                    (maxValues.z() + minValues.z()) / 2};

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
ActsExamples::MockupSectorBuilder::buildSector(
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
        det_volumes,
    const Acts::GeometryContext& gctx) {

  if (m_cfg.NumberOfSectors > maxNumberOfSectors) {
    throw std::invalid_argument("MockupSector:Number of max sectors exceeded");
    return nullptr;
  }

  // sort the detector volumes by their radial distance (from
  // innermost---->outermost)
  std::sort(det_volumes.begin(), det_volumes.end(), [](const auto& det_vol1, const auto& det_vol2) 
                                                { return det_vol1->center().y() < det_vol2->center().y(); });


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

  auto half_phi = M_PI / m_cfg.NumberOfSectors;

  if (m_cfg.NumberOfSectors == 1) {
    half_phi = (phiB - phiA) / 2;
    sector_angle = half_phi;
  }


  std::vector<float> rmins(det_volumes.size());
  std::vector<float> rmaxs(det_volumes.size());
  std::vector<float> halfZ(det_volumes.size());
  std::vector<std::unique_ptr<Acts::CylinderVolumeBounds>>cylinderVolumesBounds(det_volumes.size());

  for(std::size_t i=0; i<det_volumes.size(); i++){
    const auto &det_vol = det_volumes[i];
    rmins[i] = det_vol->center().y() - 
               det_vol->volumeBounds().values()[1] - m_cfg.toleranceOverlap;
    rmaxs[i] = std::sqrt(std::pow(det_vol->volumeBounds().values()[0],2) + 
               std::pow(det_vol->center().y() + 
               det_vol->volumeBounds().values()[1],2)) + m_cfg.toleranceOverlap; //bla bla from pythagoreum
     halfZ[i] =  det_vol->volumeBounds().values()[2];

    cylinderVolumesBounds[i] = std::make_unique<Acts::CylinderVolumeBounds>(rmins[i], rmaxs[i],
                                                   halfZ[i], sector_angle);

  }

  const Acts::Vector3 pos = {0., 0., 0.};

  // the transfom of the cylinder volume
  Acts::AngleAxis3 rotZ(M_PI / 2, Acts::Vector3(0., 0., 1));
  auto transform = Acts::Transform3(Acts::Translation3(pos));
  transform *= rotZ;

  //create a vector for the shifted surfaces of each chamber
  std::vector<std::shared_ptr<Acts::Surface>> shifted_surfaces = {};

  //creare an array of vectors that holds all the chambers of each sector

  const int array_size =  det_volumes.size();  
  std::vector<std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>> chambersOfSectors(array_size);

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
     detectorCylinderVolumesOfSector = {};

  for (int i = 0; i < m_cfg.NumberOfSectors; i++) {

    Acts::AngleAxis3 rotation(2 * i * half_phi,
                              Acts::Vector3(0., 0., 1.));

    for (std::size_t itr = 0; itr<det_volumes.size(); itr++) {

      const auto &det_vol = det_volumes[itr];

      auto shift_vol =
          rotation * Acts::Transform3(Acts::Translation3(det_vol->center()));

      for (auto& det_surf : det_vol->surfaces()) {

        auto shift_surf = Acts::Transform3::Identity() * rotation;

        // create the shifted surfaces by creating copied surface objects
        auto strawSurfaceObject = Acts::Surface::makeShared<Acts::StrawSurface>(
            det_surf->transform(Acts::GeometryContext()), 
            det_surf->bounds().values()[0], det_surf->bounds().values()[1]);

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

          chambersOfSectors[itr].push_back(detectorVolume_sec);

          shifted_surfaces.clear();

    }  // end of detector volumes

  }  // end of number of sectors


  for(std::size_t i=0; i < cylinderVolumesBounds.size(); ++i){

    detectorCylinderVolumesOfSector.push_back(Acts::Experimental::DetectorVolumeFactory::construct(
          Acts::Experimental::defaultPortalGenerator(), gctx, "cylinder_volume",
          transform, std::move(cylinderVolumesBounds[i]),
          std::vector<std::shared_ptr<Acts::Surface>>{},
          chambersOfSectors[i],
          Acts::Experimental::allPortalsAndSurfaces()));
  }


  auto cylinderVolumesBoundsOfMother = std::make_unique<Acts::CylinderVolumeBounds>(rmins.front(), rmaxs.back(),
    *std::max_element(halfZ.begin(), halfZ.end()), sector_angle);

  // creation of the mother volume
  auto detectorVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalGenerator(), gctx, "sector", transform,
      std::move(cylinderVolumesBoundsOfMother),
      std::vector<std::shared_ptr<Acts::Surface>>{}, detectorCylinderVolumesOfSector,
      Acts::Experimental::allPortalsAndSurfaces());

  return detectorVolume;
}
