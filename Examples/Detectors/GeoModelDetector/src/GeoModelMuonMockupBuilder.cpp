// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GeoModelDetector/GeoModelMuonMockupBuilder.hpp"

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

GeoModelMuonMockupBuilder::GeoModelMuonMockupBuilder(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::unique_ptr<const Acts::TrackingGeometry>
GeoModelMuonMockupBuilder::trackingGeometry(
    const Acts::GeometryContext& gctx) const {
  ConvertedVolList_t boundingBoxes = m_cfg.volumeBoxFPVs;

  // Blue print construction for the tracking geometry
  Acts::Experimental::Blueprint::Config bpCfg;
  bpCfg.envelope[Acts::AxisDirection::AxisZ] = {20_mm, 20_mm};
  bpCfg.envelope[Acts::AxisDirection::AxisR] = {5008.87000_mm, 12_mm};
  Acts::Experimental::Blueprint root{bpCfg};
  auto& cyl = root.addCylinderContainer("MuonMockupBarrelContainer",
                                        Acts::AxisDirection::AxisR);
  cyl.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
  cyl.setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

  if (boundingBoxes.empty()) {
    THROW_EXCEPTION(
        "No converted bounding boxes in the configuration - provide volumes "
        "(e.g from the GeModelDetectorObjectFactory) ");
  }

  // Add the station nodes as static cylidner nodes
  std::size_t layerId = 0;

  for (const auto& str : m_cfg.stationNames) {
    Acts::GeometryIdentifier geoIdNode =
        Acts::GeometryIdentifier().withLayer(++layerId);
    auto node = buildBarrelNode(boundingBoxes, str, *m_cfg.volumeBoundFactory,
                                geoIdNode);
    cyl.addChild(std::move(node));
  }

  auto trackingGeometry = root.construct({}, gctx, *m_logger);

  return trackingGeometry;
}

std::shared_ptr<Acts::Experimental::StaticBlueprintNode>
GeoModelMuonMockupBuilder::buildBarrelNode(
    const ConvertedVolList_t& boundingBoxes, const std::string& name,
    Acts::VolumeBoundFactory& boundFactory,
    const Acts::GeometryIdentifier& geoId) const {
  using enum Acts::TrapezoidVolumeBounds::BoundValues;

  /** Assume a station paradigm. MDT multilayers and complementary strip
   * detectors are residing under a common parent node representing a muon
   * station envelope. Group the passed boxes under by their parent */
  std::map<const GeoVPhysVol*, ConvertedVolList_t> commonStations{};
  for (const auto& box : boundingBoxes) {
    ACTS_VERBOSE("Test whether " << box.name << " contains '" << name
                                 << "' as substring");
    if (box.name.find(name) == std::string::npos) {
      continue;  // skip boxes that do not match the station name
    }
    auto parent = box.fullPhysVol->getParent().get();

    if (parent == nullptr) {
      THROW_EXCEPTION("No parent found for " << name);
    }
    commonStations[parent].push_back(box);
  }
  // Create a vector to hold the chambers
  std::vector<std::unique_ptr<Acts::TrackingVolume>> volChambers;

  if (commonStations.empty()) {
    throw std::invalid_argument("No barrel stations could be found.");
  }
  volChambers.reserve(commonStations.size());
  std::size_t stationNum = 0;
  double maxZ = std::numeric_limits<double>::lowest();
  for (const auto& [parentPhysVol, childrenTrkVols] : commonStations) {
    std::shared_ptr<Acts::Volume> parentVolume = Acts::GeoModel::convertVolume(
        Acts::GeoModel::volumePosInSpace(parentPhysVol),
        parentPhysVol->getLogVol()->getShape(), boundFactory);

    auto chamberVolume = std::make_unique<Acts::TrackingVolume>(
        *parentVolume, name + "Chamber_" + std::to_string(stationNum++));
    chamberVolume->assignGeometryId(geoId.withVolume(stationNum));

    ACTS_VERBOSE("Boundaries of the chamber volume: "
                 << chamberVolume->boundarySurfaces().size());

    std::size_t childVol = 0;
    for (const auto& child : childrenTrkVols) {
      auto trVol =
          std::make_unique<Acts::TrackingVolume>(*child.volume, child.name);
      trVol->assignGeometryId(
          geoId.withVolume(stationNum).withExtra(++childVol));

      // add the sensitives (tubes) in the constructed tracking volume
      for (const auto& surface : child.surfaces) {
        trVol->addSurface(surface);
      }

      chamberVolume->addVolume(std::move(trVol));
    }
    volChambers.push_back(std::move(chamberVolume));
    maxZ = std::max(
        maxZ, volChambers.back()->center().z() +
                  volChambers.back()->volumeBounds().values()[eHalfLengthY]);
  }

  const Acts::Vector3& cent{volChambers.front()->center()};
  double rmincyl =
      Acts::fastHypot(cent.x(), cent.y()) -
      volChambers.front()->volumeBounds().values()[eHalfLengthXnegY];
  double rmaxcyl = Acts::fastHypot(
      rmincyl +
          2 * volChambers.front()->volumeBounds().values()[eHalfLengthXnegY],
      volChambers.front()->volumeBounds().values()[eHalfLengthXposY]);
  double halfZ = maxZ;

  // Create the barrel node with the attached cylinder volume
  auto barrelNode = std::make_shared<Acts::Experimental::StaticBlueprintNode>(
      std::make_unique<Acts::TrackingVolume>(
          Acts::Transform3::Identity(),
          std::make_shared<Acts::CylinderVolumeBounds>(rmincyl, rmaxcyl, halfZ),
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
