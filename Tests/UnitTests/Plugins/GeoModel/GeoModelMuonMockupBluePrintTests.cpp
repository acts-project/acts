// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

#include <fstream>
using namespace Acts;
using namespace Acts::GeoModelReader;
using namespace Acts::Experimental;
using namespace Acts::UnitLiterals;

namespace Acts::Test {

// This tests creates the blueprint for a mockup geometry of the muon system
//  in the ATLAS experiment. The geometry is read from a SQLite database
//  and the blueprint is built using the StaticBlueprintNode class.
//  The test checks the number of children in each barrel cylinder and uses
//  static blueprint ndoes
// The following diagram shows the hierarchy of the geometry:
/*
                     +------+
                     | root |
                     +------+
                        |
         +--------------+--------------+--------------+
         |                             |              |
         v                             v              v
     +-------------+           +--------------+   +--------------+
     | InnerBarrel |           | MiddleBarrel |   | OuterBarrel  |
     +-------------+           +--------------+   +--------------+
        |                         |                  |
        v                         v                  v
     +-----------+             +-----------+       +-----------+
     | Chamber1,2 |            | Chamber1,2 |      | Chamber1,2 |
     +-----------+             +-----------+       +-----------+
        |                         |                  |
        v                         v                  v
     +-----------------+     +-----------------+   +-----------------+
     |  MultiLayer1,2  |     |  MultiLayer1,2  |   |  MultiLayer1,2  |
     +-----------------+     +-----------------+   +-----------------+
*/

constexpr std::size_t nSectors = 8;

using SensitiveSurfaces = std::vector<GeoModelSensitiveSurface>;
using BoundingBoxesVols = std::vector<
    std::pair<Acts::GeoModelDetectorObjectFactory::GeoModelBoundingBox,
              Acts::GeoModelDetectorObjectFactory::GeoModelVolumeBox>>;

// function that builds the chamber tracking volumes
std::shared_ptr<StaticBlueprintNode> buildBarrelNode(
    const BoundingBoxesVols& boundingBoxes,
    const SensitiveSurfaces& sensitiveSurfaces, const std::string& name,
    const Logger& logger) {
  // the vector with the chamber nodes and the volumes
  std::vector<std::shared_ptr<StaticBlueprintNode>> volChamberNodes;
  std::vector<std::unique_ptr<Acts::TrackingVolume>> volChambers;

  for (std::size_t sector = 0; sector < nSectors; sector++) {
    ACTS_DEBUG("Barrel name: " << name << " sector: " << sector);
    // Find the bounding box for the given chamber name and sector
    auto it_first =
        std::ranges::find_if(boundingBoxes, [&name, &sector](const auto& pair) {
          return pair.first->name().find(name) != std::string::npos &&
                 pair.first->name().find("_1_" + std::to_string(sector + 1)) !=
                     std::string::npos;
        });

    // Check if the bounding box was found
    BOOST_CHECK(it_first != boundingBoxes.end());

    auto it_second = std::find_if(
        it_first, boundingBoxes.end(), [&name, &sector](const auto& pair) {
          return pair.first->name().find(name) != std::string::npos &&
                 pair.first->name().find("_1_" + std::to_string(sector + 1)) !=
                     std::string::npos;
        });
    // Check if the second bounding box was found
    BOOST_CHECK(it_second != boundingBoxes.end());

    auto vol1 = it_first->second;
    auto vol2 = it_second->second;
    Vector3 center = (vol1->center() + vol2->center()) / 2;
    Acts::Transform3 volTrf = vol1->transform().rotation() *
                              Acts::Transform3(Acts::Translation3(center));
    auto volumeBounds = vol1->volumeBounds().values();
    double ymin = std::min(vol1->center().y(), vol2->center().y());
    double ymax = std::max(vol1->center().y(), vol2->center().y());

    double yup = ymax + volumeBounds[1];
    double ylow = ymin - volumeBounds[1];

    std::shared_ptr<Acts::CuboidVolumeBounds> chamberBounds =
        std::make_shared<Acts::CuboidVolumeBounds>(
            volumeBounds[0], (yup - ylow), volumeBounds[2]);

    std::unique_ptr<Acts::TrackingVolume> chamberVolume =
        std::make_unique<Acts::TrackingVolume>(
            volTrf, chamberBounds,
            "Chamber_" + name + "_" + std::to_string(sector + 1));

    // add the multilayer volumes
    chamberVolume->addVolume(std::make_unique<TrackingVolume>(
        *vol1, "Chamber_" + name + "_" + std::to_string(sector + 1) + "_MDT1"));
    chamberVolume->addVolume(std::make_unique<TrackingVolume>(
        *vol2, "Chamber_" + name + "_" + std::to_string(sector + 1) + "_MDT2"));

    // Now also add the RPC plane surfaces to the chamber volume
    for (const auto& surface : sensitiveSurfaces) {
      auto sname = std::get<0>(surface)->databaseEntryName();

      // Skip MDT surfaces and surfaces not in the current sector
      if (sname.find("MDT") != std::string::npos ||
          sname.find("_1_" + std::to_string(sector + 1)) == std::string::npos) {
        continue;
      }

      chamberVolume->addSurface(std::get<1>(surface));
    }

    volChambers.push_back(std::move(chamberVolume));
  }

  // create the cylinder bounds for the barrel cylinder node
  auto maxYChamberIt =
      std::max_element(volChambers.begin(), volChambers.end(),
                       [](const std::unique_ptr<Acts::TrackingVolume>& a,
                          const std::unique_ptr<Acts::TrackingVolume>& b) {
                         return a->transform().translation().y() <
                                b->transform().translation().y();
                       });

  double rmin = (*maxYChamberIt)->transform().translation().y() -
                (*maxYChamberIt)->volumeBounds().values()[1];
  double rmax = (*maxYChamberIt)->transform().translation().y() +
                (*maxYChamberIt)->volumeBounds().values()[1];
  // Create the barrel node with the attached cylinder volume
  auto barrelNode = std::make_shared<StaticBlueprintNode>(
      std::make_unique<Acts::TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(rmin, rmax, 4_m),
          name + "_Barrel"));

  // create the bluprint nodes for the chambers and add them as children to the
  // cylinder barrel node
  for (auto& chamber : volChambers) {
    auto chamberNode =
        std::make_shared<StaticBlueprintNode>(std::move(chamber));
    barrelNode->addChild(std::move(chamberNode));
  }

  return barrelNode;
}

BOOST_AUTO_TEST_SUITE(GeoModelMuonMockupBluePrintTests)

BOOST_AUTO_TEST_CASE(MockupMuonGen3Geometry) {
  // Read the geomodel file to create the mockup geometry

  const GeometryContext context;
  auto logger = Acts::getDefaultLogger("UnitTests", Acts::Logging::VERBOSE);

  Acts::GeoModelTree gmTree =
      readFromDb("/eos/user/c/cimuonsw/GeometryFiles/DimitrasMockup.db");
  std::vector<std::string> queries = {"Muon"};
  std::vector<std::unique_ptr<TrackingVolume>> chamberVolumes;

  Acts::GeoModelDetectorObjectFactory::Options options;
  options.queries = queries;

  // Create the detector object factory
  Acts::GeoModelDetectorObjectFactory::Config cfg;
  cfg.nameList = {"RPC", "Tube", "BIL", "BML", "BOL"};
  cfg.materialList = {"RPCgas", "Aluminium"};
  cfg.convertSubVolumes = true;
  cfg.convertBox = {"MDT"};

  Acts::GeoModelDetectorObjectFactory factory(cfg);
  Acts::GeoModelDetectorObjectFactory::Cache cache;
  factory.construct(cache, context, gmTree, options);
  std::cout << "box size=" << cache.boundingBoxes.size() << std::endl;

  auto innerBarrel =
      buildBarrelNode(cache.boundingBoxes, cache.sensitiveSurfaces, "BIL",
                      *logger->clone(std::nullopt, Logging::DEBUG));
  auto middleBarrel =
      buildBarrelNode(cache.boundingBoxes, cache.sensitiveSurfaces, "BML",
                      *logger->clone(std::nullopt, Logging::DEBUG));
  auto outerBarrel =
      buildBarrelNode(cache.boundingBoxes, cache.sensitiveSurfaces, "BOL",
                      *logger->clone(std::nullopt, Logging::DEBUG));

  // Check the number of children in each barrel
  BOOST_CHECK(innerBarrel->children().size() == nSectors);
  BOOST_CHECK(middleBarrel->children().size() == nSectors);
  BOOST_CHECK(outerBarrel->children().size() == nSectors);

  // Create the blueprint
  Blueprint::Config bpCfg;
  bpCfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  bpCfg.envelope[AxisDirection::AxisR] = {2_mm, 20_mm};
  Blueprint root{bpCfg};
  auto& cyl =
      root.addCylinderContainer("CylinderContainer", AxisDirection::AxisR);
  cyl.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  // Add the volumes to the blueprint
  cyl.addChild(std::move(innerBarrel));
  cyl.addChild(std::move(middleBarrel));
  cyl.addChild(std::move(outerBarrel));

  auto trackingGeometry = root.construct({}, context, *logger);
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test
