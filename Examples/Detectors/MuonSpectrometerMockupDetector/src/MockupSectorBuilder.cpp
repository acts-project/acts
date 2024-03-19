// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/MuonSpectrometerMockupDetector/MockupSectorBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/MultiWireStructureBuilder.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdaters.hpp"
#include "Acts/Plugins/Geant4/Geant4Converters.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

ActsExamples::MockupSectorBuilder::MockupSectorBuilder(
    const ActsExamples::MockupSectorBuilder::Config& config) {
  mCfg = config;
  ActsExamples::GdmlDetectorConstruction geo_gdml(mCfg.gdmlPath);
  g4World = geo_gdml.Construct();
}

std::shared_ptr<Acts::Experimental::DetectorVolume>
ActsExamples::MockupSectorBuilder::buildChamber(
    const ActsExamples::MockupSectorBuilder::ChamberConfig& chamberConfig) {
  if (g4World == nullptr) {
    throw std::invalid_argument("MockupSector: No g4World initialized");
  }

  const Acts::GeometryContext gctx;

  // Geant4Detector Config creator with the g4world from the gdml file
  auto g4WorldConfig = ActsExamples::Geant4::Geant4Detector::Config();
  g4WorldConfig.name = "Chamber";
  g4WorldConfig.g4World = g4World;

  // Get the sensitive and passive surfaces and pass to the g4World Config
  auto g4Sensitive =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          chamberConfig.SensitiveNames);
  auto g4Passive =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          chamberConfig.PassiveNames);

  auto g4SurfaceOptions = Acts::Geant4DetectorSurfaceFactory::Options();
  g4SurfaceOptions.sensitiveSurfaceSelector = g4Sensitive;
  g4SurfaceOptions.passiveSurfaceSelector = g4Passive;
  g4WorldConfig.g4SurfaceOptions = g4SurfaceOptions;

  auto g4detector = ActsExamples::Geant4::Geant4Detector();

  auto [detector, surfaces, detectorElements] =
      g4detector.constructDetector(g4WorldConfig, Acts::getDummyLogger());

  // The vector that holds the converted sensitive surfaces of the chamber
  std::vector<std::shared_ptr<Acts::Surface>> strawSurfaces = {};

  // Convert the physical volumes of the detector elements to straw surfaces
  for (auto& detectorElement : detectorElements) {
    auto context = Acts::GeometryContext();
    auto g4conv = Acts::Geant4PhysicalVolumeConverter();

    //pretranslate the detector element along z
    auto center = detectorElement->transform(gctx).translation();
    auto detElTransform = detectorElement->transform(gctx);
    detElTransform.pretranslate(Acts::Vector3(center.x(), center.y(), center.z() + mCfg.zOffset));
    //detectorElement->transform(gctx).pretranslate(Acts::Vector3(center.x(), center.y(), center.z() + mCfg.zOffset));

    g4conv.forcedType = Acts::Surface::SurfaceType::Straw;
    auto g4ConvSurf = g4conv.Geant4PhysicalVolumeConverter::surface(
        detectorElement->g4PhysicalVolume(), detElTransform);

    strawSurfaces.push_back(g4ConvSurf);
  }

  // sort the surfaces -place them in the two multilayers
  std::sort(strawSurfaces.begin(), strawSurfaces.end(),
            [&gctx](const auto& surf1, const auto& surf2) {
              if (surf1->center(gctx).x() != surf2->center(gctx).x()) {
                return surf1->center(gctx).x() < surf2->center(gctx).x();
              }
              if (surf1->center(gctx).y() != surf2->center(gctx).y()) {
                return surf1->center(gctx).y() < surf2->center(gctx).y();
              }
              return surf1->center(gctx).z() < surf2->center(gctx).z();
            });

  // split the straw surfaces for the two multilayers
  std::size_t halfSize = 0.5 * strawSurfaces.size();
  std::vector<std::shared_ptr<Acts::Surface>> strawSurfaces1(
      strawSurfaces.begin(), strawSurfaces.begin() + halfSize);
  std::vector<std::shared_ptr<Acts::Surface>> strawSurfaces2(
      strawSurfaces.begin() + halfSize, strawSurfaces.end());

  // Create the bounds of the detector volumes and the multilayers
  float radius = strawSurfaces.front()->bounds().values()[0];

  Acts::ActsScalar hx =
      strawSurfaces.front()->bounds().values()[1] + mCfg.toleranceOverlap;
  Acts::ActsScalar hy =
      0.5 * (strawSurfaces.back()->center(gctx).y() -
             strawSurfaces.front()->center(gctx).y() + 2 * radius) +
      mCfg.toleranceOverlap;
  Acts::ActsScalar hz =
      0.5 * (strawSurfaces.back()->center(gctx).z() -
             strawSurfaces.front()->center(gctx).z() + 2 * radius) +
      mCfg.toleranceOverlap;

  auto detectorVolumeBounds =
      std::make_shared<Acts::CuboidVolumeBounds>(hx, hy, hz);

  Acts::Vector3 chamber_position = {
      0.5 * (strawSurfaces.front()->center(gctx).x() +
             strawSurfaces.back()->center(gctx).x()),
      0.5 * (strawSurfaces.front()->center(gctx).y() +
             strawSurfaces.back()->center(gctx).y()),
      0.5 * (strawSurfaces.front()->center(gctx).z() + 
            strawSurfaces.back()->center(gctx).z())};

  // build the multilayers
  auto multilayer1 = buildMultiLayer(gctx, strawSurfaces1);
  auto multilayer2 = buildMultiLayer(gctx, strawSurfaces2);
  // create the detector volume for the chamber
  auto detectorVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
      chamberConfig.name,
      Acts::Transform3(Acts::Translation3(chamber_position)),
      std::move(detectorVolumeBounds), strawSurfaces,
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>{
          multilayer1, multilayer2},
      Acts::Experimental::tryAllSubVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  return detectorVolume;
}

std::shared_ptr<Acts::Experimental::DetectorVolume>
ActsExamples::MockupSectorBuilder::buildSector(
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
        detVolumes) {
  using surfacePtr = std::shared_ptr<Acts::Surface>;
  using volumePtr = std::shared_ptr<Acts::Experimental::DetectorVolume>;
  if (mCfg.NumberOfSectors > maxNumberOfSectors) {
    throw std::invalid_argument("MockupSector:Number of max sectors exceeded");
  }

  const Acts::GeometryContext gctx;

  // sort the detector volumes by their radial distance (from
  // innermost---->outermost)
  std::sort(detVolumes.begin(), detVolumes.end(),
            [](const auto& detVol1, const auto& detVol2) {
              return detVol1->center().y() < detVol2->center().y();
            });

  auto xA = detVolumes.back()->center().x() +
            detVolumes.back()->volumeBounds().values()[0];
  auto yA = detVolumes.back()->center().y() -
            detVolumes.back()->volumeBounds().values()[1];
  auto zA = detVolumes.back()->center().z();

  auto xB = detVolumes.back()->center().x() -
            detVolumes.back()->volumeBounds().values()[0];
  auto yB = detVolumes.back()->center().y() -
            detVolumes.back()->volumeBounds().values()[1];
  auto zB = detVolumes.back()->center().z();

  Acts::Vector3 pointA = {xA, yA, zA};
  Acts::Vector3 pointB = {xB, yB, zB};

  // calculate the phi angles of the vectors
  auto phiA = Acts::VectorHelpers::phi(pointA);
  auto phiB = Acts::VectorHelpers::phi(pointB);
  auto sectorAngle = M_PI;

  auto halfPhi = M_PI / mCfg.NumberOfSectors;

  if (mCfg.NumberOfSectors == 1) {
    halfPhi = (phiB - phiA) / 2;
    sectorAngle = halfPhi;
  }

  const int detVolumesSize = detVolumes.size();

  std::vector<float> rmins(detVolumesSize);
  std::vector<float> rmaxs(detVolumesSize);
  std::vector<float> halfZ(detVolumesSize);
  std::vector<std::shared_ptr<Acts::CylinderVolumeBounds>>
      cylinderVolumesBounds(detVolumesSize);

  for (int i = 0; i < detVolumesSize; i++) {
    const auto& detVol = detVolumes[i];
    rmins[i] = detVol->center().y() - detVol->volumeBounds().values()[1] -
               mCfg.toleranceOverlap;
    rmaxs[i] = std::sqrt(std::pow(detVol->volumeBounds().values()[0], 2) +
                         std::pow(detVol->center().y() +
                                      detVol->volumeBounds().values()[1],
                                  2)) +
               mCfg.toleranceOverlap;
    halfZ[i] = detVol->volumeBounds().values()[2];

    cylinderVolumesBounds[i] = std::make_shared<Acts::CylinderVolumeBounds>(
        rmins[i], rmaxs[i], halfZ[i], sectorAngle);
  }

  const Acts::Vector3 pos = {0., 0., mCfg.zOffset};

  // the transform of the cylinder volume
  Acts::AngleAxis3 rotZ(M_PI / 2, Acts::Vector3(0., 0., 1));
  auto transform = Acts::Transform3(Acts::Translation3(pos));
  transform *= rotZ;

  // create vector of surfaces
  std::vector<surfacePtr> shiftedSurfaces = {};

  // creare a vector of vectors that holds all the chambers of each sector
  std::vector<std::vector<volumePtr>> chambersOfSectors(detVolumesSize);

  std::vector<volumePtr> detectorCylinderVolumesOfSector = {};

  std::vector<volumePtr> multiLayerVolumes = {};

  // loop over the sectors
  for (int i = 0; i < mCfg.NumberOfSectors; i++) {
    Acts::AngleAxis3 rotation(2 * i * halfPhi, Acts::Vector3(0., 0., 1.));
    // loop over the layers of the chambers
    for (int itr = 0; itr < detVolumesSize; itr++) {
      const auto& detVol = detVolumes[itr];

      auto shift_vol =
          rotation * Acts::Transform3(Acts::Translation3(detVol->center()));

      for (auto& iVol : detVol->volumes()) {
        for (auto& detSurf : iVol->surfaces()) {
          auto shift_surf = Acts::Transform3::Identity() * rotation;

          // create the shifted surfaces by creating copied surface objects
          auto strawSurfaceObject =
              Acts::Surface::makeShared<Acts::StrawSurface>(
                  detSurf->transform(Acts::GeometryContext()),
                  detSurf->bounds().values()[0], detSurf->bounds().values()[1]);

          auto copiedTransformStrawSurface =
              Acts::Surface::makeShared<Acts::StrawSurface>(
                  Acts::GeometryContext(), *strawSurfaceObject, shift_surf);
          copiedTransformStrawSurface->assignGeometryId(
              Acts::GeometryIdentifier{}
                  .setExtra((int)mCfg.zOffset%255)
                  .setLayer(itr + 1)
                  .setVolume(chambersOfSectors[itr].size() + 1)
                  .setBoundary(multiLayerVolumes.size() + 1)
                  .setSensitive(shiftedSurfaces.size() + 1));

          shiftedSurfaces.push_back(copiedTransformStrawSurface);
        }  // loop over the surfaces

        auto mlVol = buildMultiLayer(gctx, shiftedSurfaces);
        mlVol->assignGeometryId(
            Acts::GeometryIdentifier{}
                .setExtra((int)mCfg.zOffset%255+1)
                .setLayer(itr + 1)
                .setVolume(chambersOfSectors[itr].size() + 1)
                .setBoundary(multiLayerVolumes.size() + 1));
        multiLayerVolumes.push_back(mlVol);
        shiftedSurfaces.clear();
      }  // loop over the multilayers

      // create the bounds of the volumes of each chamber
      auto bounds = std::make_unique<Acts::CuboidVolumeBounds>(
          detVol->volumeBounds().values()[0],
          detVol->volumeBounds().values()[1],
          detVol->volumeBounds().values()[2]);
      // create the shifted chamber
      auto detectorVolumeSec =
          Acts::Experimental::DetectorVolumeFactory::construct(
              Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
              "detectorVolumeChamber_" + std::to_string(itr) + "Sector_" +
                  std::to_string(i) +"z_" + std::to_string(mCfg.zOffset),
              shift_vol, std::move(bounds), std::vector<surfacePtr>{},
              multiLayerVolumes, Acts::Experimental::tryAllSubVolumes(),
              Acts::Experimental::tryAllPortalsAndSurfaces());

      detectorVolumeSec->assignGeometryId(
          Acts::GeometryIdentifier{}.setExtra((int)mCfg.zOffset%255).setLayer(itr + 1).setVolume(
              chambersOfSectors[itr].size() + 1));

      chambersOfSectors[itr].push_back(detectorVolumeSec);

      multiLayerVolumes.clear();

    }  // end of chambers - detector volumesdetectorVolumeChamber

  }  // end of number of sectors

  for (std::size_t i = 0; i < cylinderVolumesBounds.size(); ++i) {
    detectorCylinderVolumesOfSector.push_back(
        Acts::Experimental::DetectorVolumeFactory::construct(
            Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
            "cylinder_volume_z_" + std::to_string(mCfg.zOffset) + "_layer_" + std::to_string(i), transform,
            std::move(cylinderVolumesBounds[i]),
            std::vector<std::shared_ptr<Acts::Surface>>{}, chambersOfSectors[i],
            Acts::Experimental::tryAllSubVolumes(),
            Acts::Experimental::tryAllPortalsAndSurfaces()));
    detectorCylinderVolumesOfSector[i]->assignGeometryId(
        Acts::GeometryIdentifier{}.setExtra((int)mCfg.zOffset%255+1).setLayer(i + 1));

  }  // end of cylinder volumes

  auto cylinderVolumesBoundsOfMother =
      std::make_shared<Acts::CylinderVolumeBounds>(
          rmins.front() - mCfg.toleranceOverlap, rmaxs.back(),
          *std::max_element(halfZ.begin(), halfZ.end()), sectorAngle);

  // creation of the mother volume that captures all the volumes in phi
  auto detectorVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
      "detectorVolumeSector_z=" + std::to_string(mCfg.zOffset), transform,
      std::move(cylinderVolumesBoundsOfMother),
      std::vector<std::shared_ptr<Acts::Surface>>{},
      detectorCylinderVolumesOfSector, Acts::Experimental::tryAllSubVolumes(),
      Acts::Experimental::tryAllPortals());

   

  detectorVolume->assignGeometryId(
      Acts::GeometryIdentifier{}.setVolume((int)mCfg.zOffset%4095 + 1));

  return detectorVolume;
}

void ActsExamples::MockupSectorBuilder::drawSector(
    const std::shared_ptr<Acts::Experimental::DetectorVolume>&
        detectorVolumeSector,
    const std::string& nameObjFile) {
  Acts::ViewConfig sConfig = Acts::s_viewSensitive;

  Acts::ObjVisualization3D objSector;

  Acts::GeometryView3D::drawDetectorVolume(
      objSector, *detectorVolumeSector, Acts::GeometryContext(),
      Acts::Transform3::Identity(), sConfig);

  objSector.write(nameObjFile);
}

std::shared_ptr<Acts::Experimental::DetectorVolume>
ActsExamples::MockupSectorBuilder::buildMultiLayer(
    const Acts::GeometryContext& gctx,
    std::vector<std::shared_ptr<Acts::Surface>> surfaces) {
  // calculate the transform, bounds and the binning for the multilayer

  float radius = surfaces.front()->bounds().values()[0];
  float length = surfaces.front()->bounds().values()[1];

  Acts::Vector3 position = {0.5 * (surfaces.front()->center(gctx).x() +
                                   surfaces.back()->center(gctx).x()),
                            0.5 * (surfaces.front()->center(gctx).y() +
                                   surfaces.back()->center(gctx).y()),
                            0.5 * (surfaces.front()->center(gctx).z() +
                                   surfaces.back()->center(gctx).z())};

  float phi = acos((Acts::Vector3(0, 1, 0).dot(position)) /
                   Acts::VectorHelpers::perp(position));
  phi = (position.x() > 0) ? 2 * M_PI - phi : phi;
  Acts::AngleAxis3 rotation(phi, Acts::Vector3(0., 0., 1.));
  Acts::Transform3 mltransform =
      Acts::Transform3(Acts::Translation3(position)) * rotation;
      auto rot = mltransform.rotation();

    Acts::Vector3 rotX(rot.col(0));
    Acts::Vector3 rotY(rot.col(1));
    Acts::Vector3 rotZ(rot.col(2));

    std::cout<<rotX(0) << ", " << rotX(1)
     << ", " << rotX(2) << std::endl;
     std::cout<<rotY(0) << ", " << rotY(1)
     << ", " << rotY(2) << std::endl;
     std::cout<<rotZ(0) << ", " << rotZ(1)
     << ", " << rotZ(2) << std::endl;
     //std::cin.ignore();
  // Define the bounds in the local frame of the multilayer

  Acts::Vector3 localfront =
      mltransform.inverse() * surfaces.front()->center(gctx);
  Acts::Vector3 localback =
      mltransform.inverse() * surfaces.back()->center(gctx);

  Acts::ActsScalar hx = 0.5 * (localback.x() - localfront.x() + 2 * length);
  Acts::ActsScalar hy = 0.5 * (localback.y() - localfront.y() + 2 * radius);
  Acts::ActsScalar hz = 0.5 * (localback.z() - localfront.z() + 2 * radius);

  // apply one bin for the grid along each axis (tryAll case)
  unsigned int nBins0 = 1;
  unsigned int nBins1 = 1;

  // use the number of surfaces as the number of bins in case of binning
  if (mCfg.binning) {
    nBins0 = std::lround(hy / radius);
    nBins1 = std::lround(hz / radius);
  }

  Acts::Experimental::MultiWireStructureBuilder::Config mlCfg;

  std::unique_ptr<Acts::TrapezoidVolumeBounds> mdtBounds =
      std::make_unique<Acts::TrapezoidVolumeBounds>(hx, hx, hy, hz);

  mlCfg.name = "MultiLayer_xyz=" + std::to_string(position.x()) +
               std::to_string(position.y()) + std::to_string(position.z());
  mlCfg.mlSurfaces = surfaces;
  mlCfg.mlBounds = mdtBounds->values();
  mlCfg.transform = mltransform;
  mlCfg.mlBinning = {Acts::Experimental::ProtoBinning(
                         Acts::binY, Acts::detail::AxisBoundaryType::Bound, -hy,
                         hy, nBins0, 0u),
                     Acts::Experimental::ProtoBinning(
                         Acts::binZ, Acts::detail::AxisBoundaryType::Bound, -hz,
                         hz, nBins1, 2u)};

  Acts::Experimental::MultiWireStructureBuilder mlBuilder(mlCfg);
  return mlBuilder.construct(gctx).volumes[0];
}
