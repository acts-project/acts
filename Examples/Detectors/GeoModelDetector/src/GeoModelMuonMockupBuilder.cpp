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
#include "Acts/Utilities/Helpers.hpp"
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
  if (boundingBoxes.empty()) {
    throw std::invalid_argument(
        "GeoModelMuonMockupBuilder() -- No converted bounding boxes in the "
        "configuration - provide volumes "
        "(e.g from the GeModelDetectorObjectFactory) ");
  }

  // Blue print construction for the tracking geometry
  Acts::Experimental::Blueprint::Config bpCfg;
  bpCfg.envelope[Acts::AxisDirection::AxisZ] = {20_mm, 20_mm};
  bpCfg.envelope[Acts::AxisDirection::AxisR] = {5008.87000_mm, 12_mm};
  Acts::Experimental::Blueprint root{bpCfg};

  // Helper lambda to configure a container
  auto configureContainer = [](CylinderContainerBlueprintNode& c) {
    c.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
    c.setResizeStrategy(Acts::VolumeResizeStrategy::Gap);
  };

  // Higher level container for the muon mockup detector, stacked in Z
  auto& cyl = root.addCylinderContainer("MuonMockupContainer",
                                        Acts::AxisDirection::AxisZ);
  configureContainer(cyl);

  // First level: one container for the central cylindrical part of the MS
  // (barrel plus NSWs) and two for the Big Wheels (A and C side). These
  // containers will be stacked in Z. The Central container will have internal
  // stacking of volumes in R, while the Endcaps containers will stack station
  // volumes in Z.
  std::array<CylinderContainerBlueprintNode*,
             Acts::toUnderlying(FirstContainerIdx::nFirstContainers)>
      FirstContainers{};

  // Second level: one container for the barrel and one for the two NSWs. These
  // containers are attached to the Central container and will be stacked in R.
  // The Barrel container will have internal stacking of station volumes in R,
  // while the NSWs container will stack them in Z.
  std::array<CylinderContainerBlueprintNode*,
             Acts::toUnderlying(SecondContainerIdx::nSecondContainers)>
      SecondContainers{};

  // Helper lambda to retrieve the first-level container
  auto retrieveFirstContainer = [&cyl, &FirstContainers, &configureContainer](
                                    FirstContainerIdx containerIdx) {
    auto& container = FirstContainers[Acts::toUnderlying(containerIdx)];
    if (container == nullptr) {
      container =
          &cyl.addCylinderContainer(firstContainerIdxToString(containerIdx),
                                    (containerIdx == FirstContainerIdx::Central)
                                        ? Acts::AxisDirection::AxisR
                                        : Acts::AxisDirection::AxisZ);
      configureContainer(*container);
    }
    return container;
  };

  // Helper lambda to retrieve the second-level container
  auto retrieveSecondContainer = [&retrieveFirstContainer, &SecondContainers,
                                  &configureContainer](bool isBarrel) {
    auto& container = SecondContainers[Acts::toUnderlying(
        isBarrel ? SecondContainerIdx::Barrel : SecondContainerIdx::NSWs)];
    if (container == nullptr) {
      auto& centralContainer =
          *retrieveFirstContainer(FirstContainerIdx::Central);
      container = &centralContainer.addCylinderContainer(
          isBarrel ? "BarrelContainer" : "NSWsContainer",
          isBarrel ? Acts::AxisDirection::AxisR : Acts::AxisDirection::AxisZ);
      configureContainer(*container);
    }
    return container;
  };

  // Sorting the boxes by station
  std::ranges::sort(boundingBoxes, [this](const auto& a, const auto& b) {
    return getStationIdx(a) < getStationIdx(b);
  });

  auto it = boundingBoxes.begin();
  std::size_t layerId = 1;
  while (it != boundingBoxes.end()) {
    // Current station index
    StationIdx currentIdx = getStationIdx(*it);
    const bool isBarrel =
        (currentIdx == StationIdx ::BI || currentIdx == StationIdx ::BM ||
         currentIdx == StationIdx ::BO);

    // Find the range of boxes for the current station
    auto rangeEnd = std::find_if(it, boundingBoxes.end(),
                                 [this, currentIdx](const auto& box) {
                                   return getStationIdx(box) != currentIdx;
                                 });

    // Check we want to build this station
    const std::string station = stationIdxToString(currentIdx);
    if (!Acts::rangeContainsValue(m_cfg.stationNames, station)) {
      ACTS_DEBUG("Skipping station "
                 << station << " as not in the configured station names.");
      it = rangeEnd;
      continue;
    }
    // Build GeometryIdentifier for this station
    auto geoIdNode = Acts::GeometryIdentifier().withLayer(layerId++);

    ACTS_DEBUG("Building nodes for station " << station << " with geoId "
                                             << geoIdNode);
    auto stationNode =
        processStation(std::span(&*it, std::distance(it, rangeEnd)), station,
                       isBarrel, *m_cfg.volumeBoundFactory, geoIdNode);

    // Attach station node to the proper container
    const FirstContainerIdx firstContIdx = getFirstContainerIdx(currentIdx);
    CylinderContainerBlueprintNode* targetContainer =
        firstContIdx == FirstContainerIdx::Central
            ? retrieveSecondContainer(isBarrel)
            : retrieveFirstContainer(firstContIdx);
    targetContainer->addChild(std::move(stationNode));

    it = rangeEnd;
  }
  auto trackingGeometry = root.construct({}, gctx, *m_logger);
  if (logger().doPrint(Acts::Logging::Level::DEBUG)) {
    Acts::detail::TrackingGeometryPrintVisitor trkGeoPrinter{gctx};
    trackingGeometry->apply(trkGeoPrinter);
    ACTS_DEBUG(std::endl << trkGeoPrinter.stream().str());
  }
  return trackingGeometry;
}

GeoModelMuonMockupBuilder::NodePtr_t GeoModelMuonMockupBuilder::processStation(
    const std::span<Box_t> boundingBoxes, const std::string& station,
    const bool isBarrel, Acts::VolumeBoundFactory& boundFactory,
    const Acts::GeometryIdentifier& geoId) const {
  if (boundingBoxes.empty()) {
    ACTS_DEBUG("No chambers could be found for station" << station);
    return {};
  }
  /** Assume a station paradigm. MDT multilayers and complementary strip
   * detectors are residing under a common parent node representing a muon
   * chamber envelope. The passed boxes are grouped under their parent.
   * We create a map mapping the logical volume of each parent to its
   * StaticBlueprintNode, which contains all its children*/
  std::map<const GeoVPhysVol*, std::pair<NodePtr_t, std::size_t>>
      chamberVolumes{};
  cylBounds bounds{};
  std::size_t volNum{1};
  for (const auto& box : boundingBoxes) {
    auto parent = box.fullPhysVol->getParent().get();
    if (parent == nullptr) {
      throw std::domain_error(std::format(
          "processStation() No parent found for chamber {} in station {}",
          box.name, station));
    }

    auto it = chamberVolumes.find(parent);
    if (it == chamberVolumes.end()) {
      // We create a new chamber node for this parent
      std::shared_ptr<Acts::Volume> parentVolume =
          ActsPlugins::GeoModel::convertVolume(
              ActsPlugins::GeoModel::volumePosInSpace(parent),
              parent->getLogVol()->getShape(), boundFactory);

      auto chamberVolume = std::make_unique<Acts::TrackingVolume>(
          *parentVolume, std::format("{:}_Chamber_{:d}", station, volNum));
      chamberVolume->assignGeometryId(geoId.withVolume(volNum));

      ACTS_VERBOSE("New parent: " << chamberVolume->volumeName()
                                  << " from box: " << box.name
                                  << ", Id: " << chamberVolume->geometryId()
                                  << ", center: " << chamberVolume->center().x()
                                  << ", " << chamberVolume->center().y() << ", "
                                  << chamberVolume->center().z() << ", bounds: "
                                  << chamberVolume->volumeBounds());

      // update bounds
      if (isBarrel) {
        if (chamberVolume->volumeBounds().type() !=
            Acts::VolumeBounds::eCuboid) {
          throw std::runtime_error(std::format(
              "processStation() -- Barrel chamber {} has bound type {} instead "
              "of Cuboid",
              chamberVolume->volumeName(),
              static_cast<int>(chamberVolume->volumeBounds().type())));
        }
        updateBounds<Acts::VolumeBounds::eCuboid>(*chamberVolume, bounds);
      } else {
        if (chamberVolume->volumeBounds().type() !=
            Acts::VolumeBounds::eTrapezoid) {
          throw std::runtime_error(std::format(
              "processStation() -- Endcap chamber {} has bound type {} instead "
              "of Trapezoid",
              chamberVolume->volumeName(),
              static_cast<int>(chamberVolume->volumeBounds().type())));
        }
        updateBounds<Acts::VolumeBounds::eTrapezoid>(*chamberVolume, bounds);
      }

      it = chamberVolumes
               .emplace(parent, std::make_pair(std::make_shared<Node_t>(
                                                   std::move(chamberVolume)),
                                               volNum++))
               .first;
    }
    // We add the box to the parent node
    NodePtr_t& chamberNode = it->second.first;
    if (!chamberNode) {
      throw std::logic_error(std::format(
          "processStation() -- Found null chamber node for parent {}",
          parent->getLogVol()->getName()));
    }
    auto trVol = buildChildChamber(box, boundFactory);
    trVol->assignGeometryId(geoId.withVolume(it->second.second)
                                .withExtra(chamberNode->children().size() + 1));

    ACTS_VERBOSE("\t\t Added child: " << trVol->volumeName() << ", "
                                      << trVol->geometryId()
                                      << " to Parent: " << chamberNode->name());

    // create static blueprint node for the inner volume and add it to the
    // chamber node
    chamberNode->addChild(std::make_shared<Node_t>(std::move(trVol)));
  }
  // Create a new station node with the attached cylinder volume
  const double translationZ =
      isBarrel ? 0.0 : 0.5 * (bounds.zMax + bounds.zMin);
  const double cylHalfLenght =
      isBarrel ? bounds.zMax : 0.5 * (bounds.zMax - bounds.zMin);
  auto stationNode =
      std::make_shared<Node_t>(std::make_unique<Acts::TrackingVolume>(
          Acts::Transform3{Acts::Translation3(0., 0., translationZ)},
          boundFactory.makeBounds<Acts::CylinderVolumeBounds>(
              bounds.rMin, bounds.rMax, cylHalfLenght),
          std::format("{:}_StationVol", station)));
  ACTS_DEBUG("Created station volume: "
             << stationNode->name() << " with bounds r: " << bounds.rMin
             << " - " << bounds.rMax << ", z: " << translationZ << " pm "
             << cylHalfLenght);

  // create bluprint nodes for the chambers and add them as children to the
  // cylinder station node
  for (auto& [_, entry] : chamberVolumes) {
    stationNode->addChild(std::move(entry.first));
  }

  return stationNode;
}

std::unique_ptr<Acts::TrackingVolume>
GeoModelMuonMockupBuilder::buildChildChamber(
    const Box_t& box, Acts::VolumeBoundFactory& boundFactory) const {
  std::unique_ptr<Acts::TrackingVolume> trVol{nullptr};

  // use dedicated builder for MDT multilayers
  if (box.name.find("MDT") != std::string::npos) {
    MultiWireVolumeBuilder::Config mwCfg;
    auto vb = box.volume->volumeBoundsPtr();
    double halfY{0}, halfZ{0};
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

    mwCfg.name = box.name;
    mwCfg.mlSurfaces = box.surfaces;
    mwCfg.transform = box.volume->transform();
    auto& sb = box.surfaces.front()->bounds();
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
    trVol = std::make_unique<Acts::TrackingVolume>(*box.volume, box.name);

    // add the sensitives in the constructed tracking volume
    for (const auto& surface : box.surfaces) {
      trVol->addSurface(surface);
    }
  }
  return trVol;
}
template <Acts::VolumeBounds::BoundsType VolBounds_t>
void GeoModelMuonMockupBuilder::updateBounds(const Acts::TrackingVolume& volume,
                                             cylBounds& bounds) const {
  const Acts::Vector3& center{volume.center()};
  const double rCenter{Acts::fastHypot(center.x(), center.y())};
  const auto volBounds{volume.volumeBounds().values()};
  double rMin{0.0};
  double rMax{0.0};
  double zMin{0.0};
  double zMax{0.0};

  if constexpr (VolBounds_t == Acts::VolumeBounds::eTrapezoid) {
    using enum Acts::TrapezoidVolumeBounds::BoundValues;
    rMin = rCenter - volBounds[eHalfLengthY];
    rMax = Acts::fastHypot(rCenter + volBounds[eHalfLengthY],
                           volBounds[eHalfLengthXposY]);
    zMin = center.z() - volBounds[eHalfLengthZ];
    zMax = center.z() + volBounds[eHalfLengthZ];
  } else if constexpr (VolBounds_t == Acts::VolumeBounds::eCuboid) {
    using enum Acts::CuboidVolumeBounds::BoundValues;
    rMin = rCenter - volBounds[eHalfLengthZ];
    rMax = Acts::fastHypot(rCenter + volBounds[eHalfLengthZ],
                           volBounds[eHalfLengthX]);
    zMax = Acts::abs(center.z()) + volBounds[eHalfLengthY];
  } else {
    static_assert(VolBounds_t == Acts::VolumeBounds::eTrapezoid ||
                      VolBounds_t == Acts::VolumeBounds::eCuboid,
                  "Unsupported volume bounds type in cylBounds::update");
  }

  ACTS_VERBOSE("Computed cylindrical bounds: "
               << "r: " << rMin << ", " << rMax << " "
               << "z: " << zMin << ", " << zMax);

  bounds.rMin = std::min(bounds.rMin, rMin);
  bounds.rMax = std::max(bounds.rMax, rMax);
  bounds.zMin = std::min(bounds.zMin, zMin);
  bounds.zMax = std::max(bounds.zMax, zMax);
}

GeoModelMuonMockupBuilder::StationIdx GeoModelMuonMockupBuilder::getStationIdx(
    const Box_t& box) const {
  const auto& name = box.name;

  auto contains = [&name](std::string_view key) {
    return name.find(key) != std::string::npos;
  };
  auto isPositiveSide = [&name]() {
    // Assume stationEta is the first number following an underscore in the
    // chamber name
    std::size_t pos = name.find('_');
    while (pos != std::string::npos && pos + 1 < name.size()) {
      const char c = name[pos + 1];
      if (std::isdigit(c) != 0) {
        // we found a digit right after underscore, so positive side
        return true;
      } else if (c == '-' && pos + 2 < name.size() &&
                 std::isdigit(name[pos + 2]) != 0) {
        // we found a negative digit right after underscore, so negative side
        return false;
      }
      pos = name.find('_', pos + 1);
    }
    throw std::runtime_error("No stationEta found in name: " + name);
  };
  // Handle the inner stations
  if (contains("SmallWheel")) {
    return isPositiveSide() ? StationIdx::EAI : StationIdx::ECI;
  } else if (contains("Inner")) {
    return StationIdx::BI;
  }
  // Handle TGCs
  if (contains("TGC")) {
    return isPositiveSide() ? StationIdx::EAM : StationIdx::ECM;
  }
  // Handle MDTs and RPCs
  const bool isBarrel = contains("BMDT") || contains("RPC");
  if (contains("Middle")) {
    return isBarrel ? StationIdx::BM
                    : (isPositiveSide() ? StationIdx::EAM : StationIdx::ECM);
  } else if (contains("Outer")) {
    return isBarrel ? StationIdx::BO
                    : (isPositiveSide() ? StationIdx::EAO : StationIdx::ECO);
  } else {
    throw std::domain_error(
        "getStationIdx() -- Could not deduce station idx from volume name: " +
        name);
  }
}

GeoModelMuonMockupBuilder::FirstContainerIdx
GeoModelMuonMockupBuilder::getFirstContainerIdx(
    const StationIdx& stationIdx) const {
  if (stationIdx == StationIdx::EAM || stationIdx == StationIdx::EAO) {
    return FirstContainerIdx::BW_A;
  } else if (stationIdx == StationIdx::ECM || stationIdx == StationIdx::ECO) {
    return FirstContainerIdx::BW_C;
  } else {
    return FirstContainerIdx::Central;
  }
}

std::string GeoModelMuonMockupBuilder::stationIdxToString(
    const GeoModelMuonMockupBuilder::StationIdx idx) {
  switch (idx) {
    case StationIdx::BI:
      return "BI";
    case StationIdx::BM:
      return "BM";
    case StationIdx::BO:
      return "BO";
    case StationIdx::EAI:
      return "EAI";
    case StationIdx::EAM:
      return "EAM";
    case StationIdx::EAO:
      return "EAO";
    case StationIdx::ECI:
      return "ECI";
    case StationIdx::ECM:
      return "ECM";
    case StationIdx::ECO:
      return "ECO";
    default:
      throw std::domain_error(
          "stationIdxToString() -- Unexpected StationIdx value");
  }
}

std::string GeoModelMuonMockupBuilder::firstContainerIdxToString(
    const GeoModelMuonMockupBuilder::FirstContainerIdx idx) {
  switch (idx) {
    case FirstContainerIdx::Central:
      return "CentralContainer";
    case FirstContainerIdx::BW_A:
      return "BigWheelA_Container";
    case FirstContainerIdx::BW_C:
      return "BigWheelC_Container";
    default:
      throw std::domain_error(
          "firstContainerIdxToString() -- Unexpected FirstContainerIdx value");
  }
}
}  // namespace ActsExamples
