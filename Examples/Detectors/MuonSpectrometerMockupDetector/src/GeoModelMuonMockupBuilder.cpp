// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GeoModelMuonMockupBuilder/GeoModelMuonMockupBuilder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include "GeoModelKernel/throwExcept.h"

using namespace Acts::UnitLiterals;

namespace ActsExamples {

GeoModelMuonMockupBuilder::GeoModelMuonMockupBuilder(const Config& cfg, 
  std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::unique_ptr<const Acts::TrackingGeometry>
GeoModelMuonMockupBuilder::trackingGeometry(
    const Acts::GeometryContext& gctx) const {
  
  Acts::GeoModelDetectorObjectFactory::Config cfg;
  cfg.nameList = m_cfg.sensitivesNames;
  cfg.convertSubVolumes = true;
  cfg.convertBox = {"MDT"};
  cfg.materialList = {"Aluminium"};
  


  Acts::GeoModelDetectorObjectFactory factory(cfg);
  Acts::GeoModelDetectorObjectFactory::Cache cache;
  Acts::GeoModelDetectorObjectFactory::Options options;
  options.queries = {"Muon"};
  factory.construct(cache, gctx, m_cfg.geoModel, options);
  // Get the bounding boxes from the cache to build the geometry
  GeoModelVolumeFPVsVec boundingBoxes = cache.volumeBoxFPVs;

  // Blue print construction for the tracking geometry
  Acts::Experimental::Blueprint::Config bpCfg;
  bpCfg.envelope[Acts::AxisDirection::AxisZ] = {20_mm, 20_mm};
  bpCfg.envelope[Acts::AxisDirection::AxisR] = {2_mm, 20_mm};
  Acts::Experimental::Blueprint root{bpCfg};
  auto& cyl = root.addCylinderContainer("MuonMockupBarrelContainer",
                                        Acts::AxisDirection::AxisR);
  cyl.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
  // std::cout<<" THIS IS CRUSH"<<std::endl;
  // std::cin.ignore();
  if(boundingBoxes.empty()) {
    THROW_EXCEPTION("No bounding boxes found ");
  }

  // Add the station nodes as static cylidner nodes
  for (const auto& str : m_cfg.stationNames) {
    // std::cout<<" THIS IS CRUSH"<<std::endl;
    // std::cout<<"bpundix boxes size: "<<boundingBoxes.size()<<std::endl;
    
    auto node =
        buildBarrelNode(boundingBoxes, str, *m_cfg.volumeBoundFactory);
    cyl.addChild(std::move(node));
    // std::cout<<" THIS IS CRUSH"<<std::endl;
    // std::cin.ignore();
  }

  auto trackingGeometry = root.construct({}, gctx, *m_logger);

  return trackingGeometry;
}

std::shared_ptr<Acts::Experimental::StaticBlueprintNode>
GeoModelMuonMockupBuilder::buildBarrelNode(
    const GeoModelVolumeFPVsVec& boundingBoxes, const std::string& name,
    Acts::VolumeBoundFactory& boundFactory) const {
  /** Assume a station paradigm. MDT multilayers and complementary strip
   * detectors are residing under a common parent node representing a muon
   * station envelope. Group the passed boxes under by their parent */
    //  std::cout<<" THIS IS CRUSH"<<std::endl;
    // std::cin.ignore();
  std::map<PVConstLink, GeoModelVolumeFPVsVec> commonStations{};
 // std::cout<<" THIS IS CRUSH"<<std::endl;
  std::cout<<"bounding box size="<<boundingBoxes.size()<<std::endl;
  std::cin.ignore();
  for (const auto& box : boundingBoxes) {
    std::cout<<" THIS IS CRUSH"<<std::endl;
    //std::cin.ignore();
    auto parent = std::get<2>(box)->getParent();
    //std::cout<<" THIS IS CRUSH"<<std::endl;
    if (!parent) {
      THROW_EXCEPTION("No parent found for " << name);
    }
    commonStations[parent].push_back(box);
  }
  std::cout<<"common stations size="<<commonStations.size()<<std::endl;
  std::cin.ignore();
  // Create a vector to hold the chambers
  std::vector<std::unique_ptr<Acts::TrackingVolume>> volChambers;
  volChambers.reserve(commonStations.size());
  std::size_t stationNum = 0;
  for (const auto& [parentPhysVol, childrenTrkVols] : commonStations) {
    std::shared_ptr<Acts::Volume> parentVolume = Acts::GeoModel::convertVolume(
        parentPhysVol->getX(), parentPhysVol->getLogVol()->getShape(),
        boundFactory);
    //std::cout<<" THIS IS CRUSH"<<std::endl;

    std::unique_ptr<Acts::TrackingVolume> chamberVolume =
        std::make_unique<Acts::TrackingVolume>(
            *parentVolume, "Chamber_" + std::to_string(stationNum++));
    chamberVolume->assignGeometryId(Acts::GeometryIdentifier{}.withVolume(stationNum));
    
    std::size_t childVol = 0;
    for (const auto& child : childrenTrkVols) {
      auto& vol = std::get<0>(child);
      std::cout<<" THIS IS CRUSH in child"<<std::endl;
      //std::cin.ignore();

      std::unique_ptr<Acts::TrackingVolume> trVol =
          std::make_unique<Acts::TrackingVolume>(*vol,
                                                 std::get<1>(child)->name());
      trVol->assignGeometryId(Acts::GeometryIdentifier{}.withVolume(stationNum).withLayer(childVol++));

      // add the sensitives (tubes) in the constructed tracking volume
      auto sensitives = std::get<1>(child)->surfacePtrs();
       std::cout<<" THIS IS ANOTHER CRUSH"<<std::endl;
      //std::cin.ignore();

      for (const auto& surface : sensitives) {
        trVol->addSurface(surface);
      }

      if(chamberVolume == nullptr) {
        std::cout<<" THIS IS CRUSH for nullptr chambervol"<<std::endl;
        //std::cin.ignore();
      }

      chamberVolume->addVolume(std::move(trVol));

    }
     volChambers.push_back(std::move(chamberVolume));
  }

  std::cout<<" THIS IS CRUSH I HOPE I FOUND IT"<<std::endl;
  //std::cin.ignore();

  const Acts::Vector3& cent{volChambers.front()->center()};
  double rmincyl = Acts::fastHypot(cent.x(), cent.y()) -
                   volChambers.front()->volumeBounds().values()[0];
  double rmaxcyl = Acts::fastHypot(
      rmincyl + 2 * volChambers.front()->volumeBounds().values()[0],
      volChambers.front()->volumeBounds().values()[1]);

  // Create the barrel node with the attached cylinder volume
  auto barrelNode = std::make_shared<Acts::Experimental::StaticBlueprintNode>(
      std::make_unique<Acts::TrackingVolume>(
          Acts::Transform3::Identity(),
          std::make_shared<Acts::CylinderVolumeBounds>(rmincyl, rmaxcyl, 4_m),
          name + "_Barrel"));

  // create the bluprint nodes for the chambers and add them as children to the
  // cylinder barrel node
  for (auto& chamber : volChambers) {
    auto chamberNode =
        std::make_shared<Acts::Experimental::StaticBlueprintNode>(
            std::move(chamber));
    barrelNode->addChild(std::move(chamberNode));
  }

  return barrelNode;
}

}  // namespace ActsExamples
