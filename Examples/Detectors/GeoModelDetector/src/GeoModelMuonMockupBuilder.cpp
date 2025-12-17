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
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/MultiWireVolumeBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Geometry/detail/TrackingGeometryPrintVisitor.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <format>
#include <typeinfo>
using namespace Acts::UnitLiterals;
using namespace Acts::Experimental;

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
    throw std::invalid_argument(
        "GeoModelMuonMockupBuilder() -- No converted bounding boxes in the "
        "configuration - provide volumes "
        "(e.g from the GeModelDetectorObjectFactory) ");
  }

  // Add the station nodes as static cylidner nodes
  std::size_t layerId = 1;

  for (const auto& str : m_cfg.stationNames) {
    auto geoIdNode = Acts::GeometryIdentifier().withLayer(layerId);
    ++layerId;
    auto node = buildBarrelNode(boundingBoxes, str, *m_cfg.volumeBoundFactory,
                                geoIdNode);
    cyl.addChild(std::move(node));
  }

  auto trackingGeometry = root.construct({}, gctx, *m_logger);
  if (logger().doPrint(Acts::Logging::Level::DEBUG)) {
    Acts::detail::TrackingGeometryPrintVisitor trkGeoPrinter{gctx};
    trackingGeometry->apply(trkGeoPrinter);
    ACTS_DEBUG(std::endl << trkGeoPrinter.stream().str());
  }

  return trackingGeometry;
}

std::shared_ptr<Acts::Experimental::StaticBlueprintNode>
GeoModelMuonMockupBuilder::buildBarrelNode(
    const ConvertedVolList_t& boundingBoxes, const std::string& name,
    Acts::VolumeBoundFactory& boundFactory,
    const Acts::GeometryIdentifier& geoId) const {
  using enum Acts::CuboidVolumeBounds::BoundValues;

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
      throw std::domain_error("buildBarrelNode() No parent found for " + name);
    }
    commonStations[parent].push_back(box);
  }
  // Create a vector to hold the chambers and inner volumes
  std::vector<std::unique_ptr<Acts::TrackingVolume>> volChambers;
  std::vector<
      std::vector<std::shared_ptr<Acts::Experimental::StaticBlueprintNode>>>
      innerVolumesNodes;
  innerVolumesNodes.resize(commonStations.size());

  if (commonStations.empty()) {
    throw std::invalid_argument("No barrel stations could be found.");
  }
  volChambers.reserve(commonStations.size());
  std::size_t stationNum = 1;
  double maxZ = std::numeric_limits<double>::lowest();
  for (const auto& [parentPhysVol, childrenTrkVols] : commonStations) {
    std::shared_ptr<Acts::Volume> parentVolume =
        ActsPlugins::GeoModel::convertVolume(
            ActsPlugins::GeoModel::volumePosInSpace(parentPhysVol),
            parentPhysVol->getLogVol()->getShape(), boundFactory);

    auto chamberVolume = std::make_unique<Acts::TrackingVolume>(
        *parentVolume, std::format("{:}_Chamber_{:d}", name, stationNum));
    chamberVolume->assignGeometryId(geoId.withVolume(stationNum));

    ACTS_VERBOSE("Boundaries of the chamber volume: "
                 << chamberVolume->boundarySurfaces().size());

    std::size_t childVol = 1;
    auto chamberId = chamberVolume->geometryId();

    for (auto& child : childrenTrkVols) {
      std::unique_ptr<Acts::TrackingVolume> trVol{nullptr};

      // use dedicated builder for MDT multilayers
      if (child.name.find("MDT") != std::string::npos) {
        MultiWireVolumeBuilder::Config mwCfg;
        auto vb = child.volume->volumeBoundsPtr();
        double halfY{0};
        double halfZ{0};
        using LineBounds = Acts::LineBounds::BoundValues;

        if (vb->type() == Acts::VolumeBounds::eTrapezoid) {
          using BoundVal = Acts::TrapezoidVolumeBounds::BoundValues;

          auto tzb = std::dynamic_pointer_cast<Acts::TrapezoidVolumeBounds>(vb);
          mwCfg.bounds = boundFactory.insert(tzb);
          halfY = tzb->get(BoundVal::eHalfLengthY);
          halfZ = tzb->get(BoundVal::eHalfLengthZ);

        } else if (vb->type() == Acts::VolumeBounds::eCuboid) {
          using BoundVal = Acts::CuboidVolumeBounds::BoundValues;

          auto cbb = std::dynamic_pointer_cast<Acts::CuboidVolumeBounds>(vb);
          mwCfg.bounds = boundFactory.insert(cbb);

          halfY = cbb->get(BoundVal::eHalfLengthY);
          halfZ = cbb->get(BoundVal::eHalfLengthZ);

        } else {
          throw std::runtime_error(
              "GeoModelMuonMockupBuilder::buildBarrelNode() - Not a trapezoid "
              "or cuboid volume bounds");
        }

        mwCfg.name = child.name;
        mwCfg.mlSurfaces = child.surfaces;
        mwCfg.transform = child.volume->transform();
        auto& sb = child.surfaces.front()->bounds();
        auto lineBounds = dynamic_cast<const Acts::LineBounds*>(&sb);
        if (lineBounds == nullptr) {
          throw std::runtime_error(
              "This MDT does not have tubes, what does it have?");
        }
        double tubeR = lineBounds->get(LineBounds::eR);
        mwCfg.binning = {
            {{Acts::AxisDirection::AxisY, Acts::AxisBoundaryType::Bound, -halfY,
              halfY, static_cast<std::size_t>(std::lround(1. * halfY / tubeR))},
             2},
            {{Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, -halfZ,
              halfZ, static_cast<std::size_t>(std::lround(1. * halfZ / tubeR))},
             1}};

        MultiWireVolumeBuilder mdtBuilder{mwCfg};
        trVol = mdtBuilder.buildVolume();

      } else {
        trVol =
            std::make_unique<Acts::TrackingVolume>(*child.volume, child.name);
        trVol->assignGeometryId(chamberId.withExtra(childVol));

        // add the sensitives (tubes) in the constructed tracking volume
        for (const auto& surface : child.surfaces) {
          trVol->addSurface(surface);
        }
      }

      trVol->assignGeometryId(chamberId.withExtra(childVol));
      ++childVol;

      auto innerNode =
          std::make_shared<Acts::Experimental::StaticBlueprintNode>(
              std::move(trVol));

      innerVolumesNodes[stationNum - 1].push_back(std::move(innerNode));
    }
    volChambers.push_back(std::move(chamberVolume));
    maxZ = std::max(
        maxZ, std::abs(volChambers.back()->center().z()) +
                  volChambers.back()->volumeBounds().values()[eHalfLengthY]);
    ++stationNum;
  }

  const Acts::Vector3& cent{volChambers.front()->center()};
  double rmincyl = Acts::fastHypot(cent.x(), cent.y()) -
                   volChambers.front()->volumeBounds().values()[eHalfLengthZ];
  double rmaxcyl = Acts::fastHypot(
      rmincyl + 2 * volChambers.front()->volumeBounds().values()[eHalfLengthZ],
      volChambers.front()->volumeBounds().values()[eHalfLengthX]);
  double halfZ = maxZ;

  // Create the barrel node with the attached cylinder volume
  auto barrelNode = std::make_shared<Acts::Experimental::StaticBlueprintNode>(
      std::make_unique<Acts::TrackingVolume>(
          Acts::Transform3::Identity(),
          std::make_shared<Acts::CylinderVolumeBounds>(rmincyl, rmaxcyl, halfZ),
          std::format("{:}_Barrel", name)));

  // create the bluprint nodes for the chambers and add them as children to the
  // cylinder barrel node
  for (std::size_t chamberNum = 0; chamberNum < volChambers.size();
       ++chamberNum) {
    auto chamberNode =
        std::make_shared<Acts::Experimental::StaticBlueprintNode>(
            std::move(volChambers[chamberNum]));

    for (auto& innerVolNode : innerVolumesNodes[chamberNum]) {
      chamberNode->addChild(std::move(innerVolNode));
    }

    barrelNode->addChild(std::move(chamberNode));
  }

  return barrelNode;
}

}  // namespace ActsExamples
