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
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "GeoModelKernel/throwExcept.h"

using namespace Acts::UnitLiterals;

namespace ActsExamples {

GeoModelMuonMockupBuilder::GeoModelMuonMockupBuilder(const Config& cfg)
    : m_cfg(cfg) {
  if (m_cfg.sensitiveSurfaces.empty()) {
    THROW_EXCEPTION("No sensitive surfaces provided");
  }
}

std::unique_ptr<const Acts::TrackingGeometry>
GeoModelMuonMockupBuilder::trackingGeometry(
    const Acts::GeometryContext& gctx) const {
  // Blue print construction for the tracking geometry
  Acts::Experimental::Blueprint::Config bpCfg;
  bpCfg.envelope[Acts::AxisDirection::AxisZ] = {20_mm, 20_mm};
  bpCfg.envelope[Acts::AxisDirection::AxisR] = {2_mm, 20_mm};
  Acts::Experimental::Blueprint root{bpCfg};
  auto& cyl = root.addCylinderContainer("MuonMockupBarrelContainer",
                                        Acts::AxisDirection::AxisR);
  cyl.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

  // Add the station nodes as static cylidner nodes
  for (const auto& str : m_cfg.stationNames) {
    auto node =
        buildBarrelNode(m_cfg.boundingBoxes, m_cfg.sensitiveSurfaces, str);
    cyl.addChild(std::move(node));
  }

  auto trackingGeometry = root.construct({}, gctx);

  return trackingGeometry;
}

std::shared_ptr<Acts::Experimental::StaticBlueprintNode>
GeoModelMuonMockupBuilder::buildBarrelNode(
    const GeoModelVolumeFPVsVec& boundingBoxes,
    const SensitiveSurfaces& sensitiveSurfaces, const std::string& name) const {
  // Create a vector to hold the chambers
  std::vector<std::unique_ptr<Acts::TrackingVolume>> volChambers;
  volChambers.reserve(m_cfg.nSectors);
  for (std::size_t sector = 0; sector < m_cfg.nSectors; sector++) {
    // ACTS_DEBUG("Barrel name: " << name << " sector: " << sector);
    // Find the bounding box for the given chamber name and sector
    auto it_first =
        std::ranges::find_if(boundingBoxes, [&name, &sector](const auto& box) {
          auto volName = std::get<1>(box)->name();
          return volName.find(name) != std::string::npos &&
                 volName.find("MDT") != std::string::npos &&
                 volName.ends_with(std::to_string(sector + 1) + "_1");
        });

    auto it_second =
        std::ranges::find_if(boundingBoxes, [&name, &sector](const auto& box) {
          auto volName = std::get<1>(box)->name();
          return volName.find(name) != std::string::npos &&
                 volName.find("MDT") != std::string::npos &&
                 volName.ends_with(std::to_string(sector + 1) + "_2");
        });

    if (it_first == boundingBoxes.end() || it_second == boundingBoxes.end()) {
      THROW_EXCEPTION("No bounding box found ");
      continue;
    }

    // print the names of the bounding boxes

    // Construct the chamber from the multilayer volumes
    Acts::Volume vol1 = std::get<0>(*it_first);
    Acts::Volume vol2 = std::get<0>(*it_second);

    // construct the tracking volumes from the multilayer volumes
    auto trVol1 = std::make_unique<Acts::TrackingVolume>(
        vol1,
        "MultiLayer_" + name + "_" + std::to_string(sector + 1) + "_MDT1");

    auto trVol2 = std::make_unique<Acts::TrackingVolume>(
        vol2,
        "MultiLayer_" + name + "_" + std::to_string(sector + 1) + "_MDT2");

    auto parent = std::get<2>(*it_first)->getParent();

    if (!parent) {
      THROW_EXCEPTION("No parent found for "<<name <<" sector: " <<sector);
    }

    // COnstruct the parent as the chamber tracking volume

    Acts::Volume parentVolume = Acts::GeoModel::convertVolume(
        parent->getX(), *parent->getLogVol()->getShape());

    std::unique_ptr<Acts::TrackingVolume> chamberVolume =
        std::make_unique<Acts::TrackingVolume>(
            parentVolume, "Chamber_" + name + "_" + std::to_string(sector + 1));

    // add the sensitive surfaces to the chamber volume (tubes and RPCs)
    for (const auto& surface : sensitiveSurfaces) {
      auto sname = std::get<0>(surface)->databaseEntryName();

      // Skip MDT surfaces and surfaces not in the current sector
      if (sname.find("MDT") != std::string::npos ||
          sname.find("RPC11_1_" + std::to_string(sector + 1)) ==
              std::string::npos ||
          sname.find(name) == std::string::npos) {
        continue;
      }
      chamberVolume->addSurface(std::get<1>(surface));
    }

    // get the MDT tubes from the multilayer detector volumes
    auto detVol1 = std::get<1>(*it_first);
    auto detVol2 = std::get<1>(*it_second);

    for (const auto& surf : detVol1->surfacePtrs()) {
      trVol1->addSurface(surf);
    }
    for (const auto& surf : detVol2->surfacePtrs()) {
      trVol2->addSurface(surf);
    }

    // Add the multilayer volumes to the chamber volume
    chamberVolume->addVolume(std::move(trVol1));
    chamberVolume->addVolume(std::move(trVol2));

    volChambers.push_back(std::move(chamberVolume));
  }

  // create the cylinder bounds for the barrel cylinder node

  double rmincyl = volChambers.front()->center().perp() -
                   volChambers.front()->volumeBounds().values()[0];
  double rmaxcyl =
      std::hypot(rmincyl + 2 * volChambers.front()->volumeBounds().values()[0],
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