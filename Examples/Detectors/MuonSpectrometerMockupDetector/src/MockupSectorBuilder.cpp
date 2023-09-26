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
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/MultiWireStructureBuilder.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/interface/IExternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Plugins/Geant4/Geant4Converters.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"
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
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

/// @brief  Internal Builder for volumes with internal volumes
 class InternalVolumesBuilder : public Acts::Experimental::IInternalStructureBuilder{
 public:

  /// Constructor
  ///
  ///@param volumes The internal volumes of the internal structure
  InternalVolumesBuilder(std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> volumes)
  : Acts::Experimental::IInternalStructureBuilder(),
  m_volumes(std::move(volumes)) {}

  Acts::Experimental::InternalStructure construct(
    [[maybe_unused]] const Acts::GeometryContext& gctx) const final {

    return {{}, m_volumes, Acts::Experimental::tryAllPortals(), Acts::Experimental::tryAllSubVolumes()};
  }

  private:

    ///  internal volumes
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> m_volumes;

};

/// @brief External Builder for volumes with internal volumes
template<typename bounds_type>
class ExternalVolumesBuilder : public Acts::Experimental::IExternalStructureBuilder{
public:

  /// Constructor
  ///
  ///@param transform The transform for the volume
  ///@param bounds The cuboid volume bounds
  ExternalVolumesBuilder(const Acts::Transform3& transform, const bounds_type& bounds)
  : Acts::Experimental::IExternalStructureBuilder(),
  m_transform(transform),
  m_bounds(std::move(bounds)) {}

  Acts::Experimental::ExternalStructure construct( 
    [[maybe_unused]] const Acts::GeometryContext& gctx) const final{

    return{ m_transform, std::make_unique<bounds_type>(m_bounds), 
    Acts::Experimental::defaultPortalAndSubPortalGenerator()};
  }

private:

  Acts::Transform3 m_transform;

  bounds_type m_bounds;


};


ActsExamples::MockupSectorBuilder::MockupSectorBuilder(
    const ActsExamples::MockupSectorBuilder::Config& config, std::unique_ptr<const Acts::Logger> logger) {
  mCfg = config;
  mLogger = std::move(logger);
  if(mCfg.gdmlPath != ""){
  ActsExamples::GdmlDetectorConstruction geo_gdml(mCfg.gdmlPath);
  g4World = geo_gdml.Construct();
}

}

 std::shared_ptr<Acts::Experimental::DetectorComponent> 
 ActsExamples::MockupSectorBuilder::buildMultiLayer(
     ActsExamples::MockupSectorBuilder::MultiLayerConfig& mlConfig){

  if (g4World == nullptr ) {
    throw std::invalid_argument("MockupSector: No g4World initialized");
    return nullptr;
  }

   const Acts::GeometryContext gctx;
     //The vector that holds the surfaces of the chamber
  std::vector<std::shared_ptr<Acts::Surface>> strawSurfaces = {};

    // Geant4Detector Config creator with the g4world from the gdml file
  auto g4WorldConfig = ActsExamples::Geant4::Geant4Detector::Config();
  g4WorldConfig.name = mlConfig.name;
  g4WorldConfig.g4World = g4World;

  // Get the sensitive and passive surfaces and pass to the g4World Config
  auto g4Sensitive =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          mlConfig.SensitiveNames);
  auto g4Passive =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          mlConfig.PassiveNames);

  auto g4SurfaceOptions = Acts::Geant4DetectorSurfaceFactory::Options();
  g4SurfaceOptions.sensitiveSurfaceSelector = g4Sensitive;
  g4SurfaceOptions.passiveSurfaceSelector = g4Passive;
  g4WorldConfig.g4SurfaceOptions = g4SurfaceOptions;

  auto g4detector = ActsExamples::Geant4::Geant4Detector();

  auto [detector, surfaces, detectorElements] =
      g4detector.constructDetector(g4WorldConfig, Acts::getDummyLogger());

  std::array<std::pair<float, float>, 3> min_max;
  std::fill(min_max.begin(), min_max.end(),
            std::make_pair<float, float>(std::numeric_limits<float>::max(),
                                         -std::numeric_limits<float>::max()));

  // Convert the physical volumes of the detector elements to straw surfaces
  for (const auto& detectorElement : detectorElements) {
   
    auto g4conv = Acts::Geant4PhysicalVolumeConverter();

    g4conv.forcedType = Acts::Surface::SurfaceType::Straw;
    auto g4ConvSurf = g4conv.Geant4PhysicalVolumeConverter::surface(
        detectorElement->g4PhysicalVolume(),
        detectorElement->transform(gctx));

    strawSurfaces.push_back(g4ConvSurf);

     min_max[0].first =
        std::min(min_max[0].first, (float)g4ConvSurf->center(gctx).x());
    min_max[0].second =
        std::max(min_max[0].second, (float)g4ConvSurf->center(gctx).x());

    min_max[1].first =
        std::min(min_max[1].first, (float)g4ConvSurf->center(gctx).y());
    min_max[1].second =
        std::max(min_max[1].second, (float)g4ConvSurf->center(gctx).y());

    min_max[2].first =
        std::min(min_max[2].first, (float)g4ConvSurf->center(gctx).z());
    min_max[2].second =
        std::max(min_max[2].second, (float)g4ConvSurf->center(gctx).z());

  }
  float radius = strawSurfaces.front()->bounds().values()[0];

  Acts::Vector3 minValues = {min_max[0].first, min_max[1].first,
                             min_max[2].first};
  Acts::Vector3 maxValues = {min_max[0].second, min_max[1].second,
                             min_max[2].second};

  Acts::ActsScalar hx =
      strawSurfaces.front()->bounds().values()[1] + mCfg.toleranceOverlap;
  Acts::ActsScalar hy =
      0.5 * ((maxValues.y() + radius) - (minValues.y() - radius)) +
      mCfg.toleranceOverlap;
  Acts::ActsScalar hz =
      0.5 * ((maxValues.z() + radius) - (minValues.z() - radius)) +
      mCfg.toleranceOverlap;

  auto ypos = 0.5 * (minValues.y() + maxValues.y());


  const Acts::Vector3 pos = {0., ypos, 0.};

  // Use the Multiwire Structure Builder to build the multilayer component
  Acts::Experimental::MultiWireStructureBuilder::Config mlwCfg;
  mlwCfg.name = mlConfig.name;
  mlwCfg.mlSurfaces = strawSurfaces;
  mlwCfg.mlBounds = {hx, hy, hz};
  mlwCfg.transform = Acts::Transform3(Acts::Translation3(pos));
  mlwCfg.mlBinning = defineBinning(std::make_pair(pos.z()-hz,pos.z()+hz), std::make_pair(pos.y()-hy, pos.y()+hy), radius);

  auto mlwBuilder = std::make_shared<Acts::Experimental::MultiWireStructureBuilder>(mlwCfg);

  auto multiLayer_component = std::make_shared<Acts::Experimental::DetectorComponent>(mlwBuilder->construct(gctx));

  return multiLayer_component;
}

std::shared_ptr<Acts::Experimental::DetectorComponent>
ActsExamples::MockupSectorBuilder::buildChamber(
    const ActsExamples::MockupSectorBuilder::ChamberConfig& chamberConfig) {

  const Acts::GeometryContext gctx;

  auto volumes = chamberConfig.internalVolumes;

  //sorting the internal volumes 
  std::sort(volumes.begin(), volumes.end(),
    [](const auto& vol1, const auto& vol2){
      return vol1->center().y() < vol2->center().y();
    });


  //Calculate the bounds of the chamber from the inernal volumes
  Acts::ActsScalar hy = 0.5*((volumes.back()->center().y() + 
    volumes.back()->volumeBounds().values()[1]) - 
  (volumes.front()->center().y() - volumes.front()->volumeBounds().values()[1]) + mCfg.toleranceOverlap);

  Acts::ActsScalar hx = volumes.front()->volumeBounds().values()[0] + mCfg.toleranceOverlap;
  Acts::ActsScalar hz = volumes.front()->volumeBounds().values()[2] + mCfg.toleranceOverlap;

  auto midx = 0.5*(volumes.back()->center().x() + volumes.front()->center().x());
  auto midy = 0.5*(volumes.back()->center().y() + volumes.front()->center().y());
  auto midz = 0.5*(volumes.back()->center().z() + volumes.front()->center().z());

  //Set the externals builder

  auto transform = Acts::Transform3(Acts::Translation3(Acts::Vector3(midx,midy,midz)));
  Acts::CuboidVolumeBounds bounds(hx,hy,hz);

  auto ebuilder = std::make_shared<ExternalVolumesBuilder<Acts::CuboidVolumeBounds>>(
    transform, bounds);

  //Set the internals builder

  auto ibuilder = std::make_shared<InternalVolumesBuilder>(volumes);

  //Construct chamber detector volume

  Acts::Experimental::DetectorVolumeBuilder::Config dvConfig;
  dvConfig.auxiliary = "Construct Detector Volume for the chamber";
  dvConfig.name = chamberConfig.name;
  dvConfig.internalsBuilder = ibuilder;
  dvConfig.externalsBuilder = ebuilder;

   auto dvBuilder = std::make_shared<Acts::Experimental::DetectorVolumeBuilder>(
      dvConfig,
      Acts::getDefaultLogger("DetectorVolumeBuilder", Acts::Logging::VERBOSE));

  auto chamber_component = std::make_shared<Acts::Experimental::DetectorComponent>(dvBuilder->construct(gctx));

  return chamber_component;

}

std::shared_ptr<const Acts::Experimental::DetectorComponent>
ActsExamples::MockupSectorBuilder::buildSector(
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
        detVolumes) {
  if (mCfg.NumberOfSectors > maxNumberOfSectors) {
    throw std::invalid_argument("MockupSector:Number of max sectors exceeded");
    return nullptr;
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

  const Acts::Vector3 pos = {0., 0., 0.};

  // the transform of the mother cylinder volume
  Acts::AngleAxis3 rotZ(M_PI / 2, Acts::Vector3(0., 0., 1));
  auto transform = Acts::Transform3(Acts::Translation3(pos));
  transform *= rotZ;

  // create a vector for the shifted surfaces of each chamber
  std::vector<std::shared_ptr<Acts::Surface>> shiftedSurfaces = {};
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> shiftedMultiLayers = {};

  // creare a vector of vectors that holds all the chambers of each sector
  std::vector<std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>>
      chambersOfSectors(detVolumesSize);

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detectorCylinderVolumesOfSector = {};

  for (int i = 0; i < mCfg.NumberOfSectors; i++) {
    Acts::AngleAxis3 rotation(2 * i * halfPhi, Acts::Vector3(0., 0., 1.));

    for (int itr = 0; itr < detVolumesSize; itr++) {
      const auto& detVol = detVolumes[itr];

      auto shift_vol =
          rotation * Acts::Transform3(Acts::Translation3(detVol->center()));

          //loop over inernal volumes (e.g multilayers)
          auto iVolumes = detVol->volumes();

        for(std::size_t iv=0; iv < iVolumes.size(); iv++){

            const auto& ivol = iVolumes[iv];

            auto shift_ivol = 
            rotation * Acts::Transform3(Acts::Translation3(ivol->center()));

        for (auto& detSurf : ivol->surfaces()) {
         auto shift_surf = Acts::Transform3::Identity() * rotation;

          // create the shifted surfaces by creating copied surface objects
          auto strawSurfaceObject = Acts::Surface::makeShared<Acts::StrawSurface>(
            detSurf->transform(Acts::GeometryContext()),
            detSurf->bounds().values()[0], detSurf->bounds().values()[1]);

          auto copiedTransformStrawSurface =
            Acts::Surface::makeShared<Acts::StrawSurface>(
                Acts::GeometryContext(), *strawSurfaceObject, shift_surf);

          shiftedSurfaces.push_back(copiedTransformStrawSurface);
      }// end of surfaces

      auto radius = shiftedSurfaces.front()->bounds().values()[0];

      //Build the shifted multi layers with the multi wire structure builder
      Acts::Experimental::MultiWireStructureBuilder::Config mlwCfg;
      mlwCfg.name = "Shifted_MultiLayer" + std::to_string(iv) + 
      "Chamber" + std::to_string(itr) + "Sector" +std::to_string(i);
      mlwCfg.mlSurfaces = shiftedSurfaces;
      mlwCfg.transform = shift_ivol;
      mlwCfg.mlBounds = {ivol->volumeBounds().values()[0], 
      ivol->volumeBounds().values()[1],
      ivol->volumeBounds().values()[2]};
      mlwCfg.mlBinning = defineBinning(std::make_pair(shift_ivol.translation().z() - ivol->volumeBounds().values()[2], shift_ivol.translation().z() + ivol->volumeBounds().values()[2]),
        std::make_pair(shift_ivol.translation().y() - ivol->volumeBounds().values()[1], shift_ivol.translation().y() + ivol->volumeBounds().values()[1]), radius);

      auto mlwBuilder = std::make_shared<Acts::Experimental::MultiWireStructureBuilder>(mlwCfg);

      auto multiLayer_component = std::make_shared<Acts::Experimental::DetectorComponent>(mlwBuilder->construct(gctx));

          shiftedMultiLayers.push_back(multiLayer_component->volumes.front());
          shiftedSurfaces.clear();
        } // end of inner volumes - multilayers

      // create the bounds of the shifted chamber
    auto bounds = std::make_unique<Acts::CuboidVolumeBounds>(
          detVol->volumeBounds().values()[0],
          detVol->volumeBounds().values()[1],
          detVol->volumeBounds().values()[2]);
      // create the shifted chamber
      auto detectorVolumeSec =
          Acts::Experimental::DetectorVolumeFactory::construct(
              Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
              "detectorVolumeChamber_" + std::to_string(itr)+"Sector_"+std::to_string(i), shift_vol,
              std::move(bounds), 
              std::vector<std::shared_ptr<Acts::Surface>>{},
              shiftedMultiLayers,
              Acts::Experimental::tryAllSubVolumes(),
              Acts::Experimental::tryAllPortals());

      chambersOfSectors[itr].push_back(detectorVolumeSec);
      shiftedMultiLayers.clear();


    }  // end of detector volumes

  }  // end of number of sectors

  for (std::size_t i = 0; i < cylinderVolumesBounds.size(); ++i) {
    detectorCylinderVolumesOfSector.push_back(
        Acts::Experimental::DetectorVolumeFactory::construct(
            Acts::Experimental::defaultPortalAndSubPortalGenerator(), gctx,
            "cylinder_volume_" + std::to_string(i), transform,
            std::move(cylinderVolumesBounds[i]),
            std::vector<std::shared_ptr<Acts::Surface>>{}, chambersOfSectors[i],
            Acts::Experimental::tryAllSubVolumes(),
            Acts::Experimental::tryAllPortals()));
  }  // end of cylinder volumes

  //Internal Builder
  auto ibuilder = std::make_shared<InternalVolumesBuilder>(detectorCylinderVolumesOfSector);

  //External Builder
  Acts::CylinderVolumeBounds cylinderVolumesBoundsOfMother(rmins.front(), rmaxs.back(),
          *std::max_element(halfZ.begin(), halfZ.end()), sectorAngle);

  auto ebuilder = std::make_shared<ExternalVolumesBuilder<Acts::CylinderVolumeBounds>>(transform,cylinderVolumesBoundsOfMother);

  Acts::Experimental::DetectorVolumeBuilder::Config dvCfg;
  dvCfg.auxiliary = "Construct Detector Volume for the Sector";
  dvCfg.name = "SectorVolume";
  dvCfg.externalsBuilder = ebuilder;
  dvCfg.internalsBuilder = ibuilder;
  //dvCfg.addInternalsToRoot = true;

  auto dvBuilder = std::make_shared<Acts::Experimental::DetectorVolumeBuilder>(dvCfg, Acts::getDefaultLogger("DetectorVolumeBuilder", Acts::Logging::VERBOSE));


  auto sector_component = std::make_shared<Acts::Experimental::DetectorComponent>(dvBuilder->construct(gctx));

  return sector_component;
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

std::vector<Acts::Experimental::ProtoBinning> ActsExamples::MockupSectorBuilder::defineBinning(std::pair<double,double> binEdgesZ, std::pair<double,double> binEdgesY, 
  float radius){

   if(mCfg.robustMode){

    ACTS_VERBOSE("Building Multi Layer with robust mode - 1 bin along y and z axis without bin expansion");
    return { Acts::Experimental::ProtoBinning(Acts::binZ, Acts::detail::AxisBoundaryType::Bound,
                                          binEdgesZ.first, binEdgesZ.second, 1),
    Acts::Experimental::ProtoBinning(Acts::binY, Acts::detail::AxisBoundaryType::Bound,
                                          binEdgesY.first, binEdgesY.second, 1)};
     
  }else{ 

      ACTS_VERBOSE("Building Multi Layer - 1 surface per bin along y and z axis with 1 bin expansion");

      auto binWidth = 2 *radius;
      auto nBinsY = (binEdgesY.second - binEdgesY.first)/binWidth;
      auto nBinsZ = (binEdgesZ.second - binEdgesZ.first)/binWidth;
      return {
    Acts::Experimental::ProtoBinning(Acts::binZ, Acts::detail::AxisBoundaryType::Bound,
      binEdgesZ.first, binEdgesZ.second, nBinsZ, 1u),
    Acts::Experimental::ProtoBinning(Acts::binY, Acts::detail::AxisBoundaryType::Bound,
      binEdgesY.first, binEdgesY.second, nBinsY, 0u)};

    }

}
