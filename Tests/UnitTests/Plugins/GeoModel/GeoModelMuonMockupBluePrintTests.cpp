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
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
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
     +--------------+             +-------------+       +-------------+
     | Chamber1,2.. |            | Chamber1,2.. |      | Chamber1,2.. |
     +--------------+             +-------------+       +-------------+
        |                         |                  |
        v                         v                  v
     +-----------------+     +-----------------+   +-----------------+
     |  MultiLayer1,2  |     |  MultiLayer1,2  |   |  MultiLayer1,2  |
     +-----------------+     +-----------------+   +-----------------+
*/

constexpr std::size_t nSectors = 8;
const GeometryContext gctx;

using SensitiveSurfaces = std::vector<GeoModelSensitiveSurface>;
using BoundingBoxesVols = std::vector<
    std::pair<Acts::GeoModelDetectorObjectFactory::GeoModelBoundingBox,
              Acts::GeoModelDetectorObjectFactory::GeoModelVolumeBox>>;

// function that builds the BluePrint Nodes for the barrel - cylinders
std::shared_ptr<StaticBlueprintNode> buildBarrelNode(
    const BoundingBoxesVols& boundingBoxes,
    const SensitiveSurfaces& sensitiveSurfaces, const std::string& name,
    const GeometryContext& context, const Logger& logger) {
  // the vector with the chamber nodes and the volumes
  std::vector<std::shared_ptr<StaticBlueprintNode>> volChamberNodes;
  std::vector<std::unique_ptr<Acts::TrackingVolume>> volChambers;

  for (std::size_t sector = 0; sector < nSectors; sector++) {
    ACTS_DEBUG("Barrel name: " << name << " sector: " << sector);
    // Find the bounding box for the given chamber name and sector
    auto it_first = std::ranges::find_if(boundingBoxes, [&name, &sector](
                                                            const auto& pair) {
      return pair.first->name().find(name) != std::string::npos &&
             pair.first->name().find("MDT03_1_" + std::to_string(sector + 1)) !=
                 std::string::npos;
    });
    for (const auto& box : boundingBoxes) {
      std::cout << "Box name: " << box.first->name() << std::endl;
    }

    // Check if the bounding box was found
    BOOST_CHECK(it_first != boundingBoxes.end());

    auto it_second = std::ranges::find_if(boundingBoxes, [&name, &sector](
                                                             const auto& pair) {
      return pair.first->name().find(name) != std::string::npos &&
             pair.first->name().find("MDT04_1_" + std::to_string(sector + 1)) !=
                 std::string::npos;
    });
    // Check if the second bounding box was found
    BOOST_CHECK(it_second != boundingBoxes.end());
    // print the names of the bounding boxes

    // Construct the chamber from the multilayer volumes
    auto vol1 = it_first->second;
    auto vol2 = it_second->second;
    // print the center of the volumes

    double r1 = std::hypot(vol1->center().x(), vol1->center().y());
    double r2 = std::hypot(vol2->center().x(), vol2->center().y());
    double rmax = std::max(r1, r2);
    double rmin = std::min(r1, r2);

    auto volumeBounds = vol1->volumeBounds().values();

    // calculate the transform and the bounds of the wrapping chamber for the
    // two multilayers
    double rOuter = rmax + volumeBounds[1];
    double rInner = rmin - volumeBounds[1];
    double rCenter = (rOuter + rInner) / 2;
    double centerZ = vol1->center().z();
    double phiAngle = std::atan2(vol1->center().y(), vol1->center().x());
    double centerX = rCenter * std::cos(phiAngle);
    double centerY = rCenter * std::sin(phiAngle);
    Acts::Vector3 center(centerX, centerY, centerZ);
    Acts::Transform3 volTrf = Acts::Transform3(Acts::Translation3(center));
    volTrf.linear() = vol1->transform().rotation();

    // Now also add the RPC plane surfaces to the chamber volume
    std::vector<std::shared_ptr<Acts::Surface>> rpcSurfaces;
    for (const auto& surface : sensitiveSurfaces) {
      auto sname = std::get<0>(surface)->databaseEntryName();

      // Skip MDT surfaces and surfaces not in the current sector
      if (sname.find("MDT") != std::string::npos ||
          sname.find("RPC11_1_" + std::to_string(sector + 1)) ==
              std::string::npos ||
          sname.find(name) == std::string::npos) {
        continue;
      }
      rpcSurfaces.push_back(std::get<1>(surface));
    }
    // sort the rpc surfaces by their radial distance
    std::sort(rpcSurfaces.begin(), rpcSurfaces.end(),
              [context](const auto& a, const auto& b) {
                return std::hypot(a->center(gctx).x(), a->center(gctx).y()) <
                       std::hypot(b->center(gctx).x(), b->center(gctx).y());
              });

    double minRpc = std::hypot(rpcSurfaces.front()->center(gctx).x(),
                               rpcSurfaces.front()->center(gctx).y());
    double maxRpc = std::hypot(rpcSurfaces.back()->center(gctx).x(),
                               rpcSurfaces.back()->center(gctx).y());

    std::shared_ptr<Acts::TrapezoidVolumeBounds> chamberBounds =
        std::make_shared<Acts::TrapezoidVolumeBounds>(
            0.5 * (maxRpc - minRpc), 0.5 * (maxRpc - minRpc), volumeBounds[2],
            volumeBounds[3]);

    std::unique_ptr<Acts::TrackingVolume> chamberVolume =
        std::make_unique<Acts::TrackingVolume>(
            volTrf, chamberBounds,
            "Chamber_" + name + "_" + std::to_string(sector + 1));

    // add the multilayer volumes as tracking volumes
    // fill the multilayer volumes with the tubes of the mdts
    std::unique_ptr<Acts::TrackingVolume> trVol1 =
        std::make_unique<Acts::TrackingVolume>(
            vol1->transform(), vol1->volumeBoundsPtr(),
            "MultiLayer_" + name + "_" + std::to_string(sector + 1) + "_MDT1");
    std::unique_ptr<Acts::TrackingVolume> trVol2 =
        std::make_unique<Acts::TrackingVolume>(
            vol2->transform(), vol2->volumeBoundsPtr(),
            "MultiLayer_" + name + "_" + std::to_string(sector + 1) + "_MDT2");

    for (const auto& surf : it_first->first->surfacePtrs()) {
      std::cout << "Adding MDT surface: " << surf->name() << std::endl;
      trVol1->addSurface(surf);
    }
    for (const auto& surf : it_second->first->surfacePtrs()) {
      std::cout << "Adding MDT surface: " << surf->name() << std::endl;
      trVol2->addSurface(surf);
    }

    chamberVolume->addVolume(std::move(trVol1));
    chamberVolume->addVolume(std::move(trVol2));

    for (const auto& surface : rpcSurfaces) {
      std::cout << "Adding RPC surface: " << surface->name() << std::endl;
      chamberVolume->addSurface(surface);
    }

    volChambers.push_back(std::move(chamberVolume));
  }

  // create the cylinder bounds for the barrel cylinder node

  double rmincyl = std::hypot(volChambers.front()->center().x(),
                              volChambers.front()->center().y()) -
                   volChambers.front()->volumeBounds().values()[1];
  double rmaxcyl =
      std::hypot(rmincyl + 2 * volChambers.front()->volumeBounds().values()[0],
                 volChambers.front()->volumeBounds().values()[3]);

  // Create the barrel node with the attached cylinder volume
  auto barrelNode = std::make_shared<StaticBlueprintNode>(
      std::make_unique<Acts::TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(rmincyl, rmaxcyl, 4_m),
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
  factory.construct(cache, gctx, gmTree, options);
  std::cout << "box size=" << cache.boundingBoxes.size() << std::endl;

  auto innerBarrel =
      buildBarrelNode(cache.boundingBoxes, cache.sensitiveSurfaces, "BIL", gctx,
                      *logger->clone(std::nullopt, Logging::DEBUG));
  auto middleBarrel =
      buildBarrelNode(cache.boundingBoxes, cache.sensitiveSurfaces, "BML", gctx,
                      *logger->clone(std::nullopt, Logging::DEBUG));
  auto outerBarrel =
      buildBarrelNode(cache.boundingBoxes, cache.sensitiveSurfaces, "BOL", gctx,
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

  BOOST_CHECK(cyl.children().size() == 3);

  auto trackingGeometry = root.construct({}, gctx, *logger);

  std::size_t nVolumes = 0;
  trackingGeometry->visitVolumes([&](const TrackingVolume* volume) {
    std::string name = volume->volumeName();
    if (name.find("Chamber") != std::string::npos) {
      nVolumes++;
    }
  });

  std::size_t nMultiLayers = 0;
  trackingGeometry->visitVolumes([&](const TrackingVolume* volume) {
    std::string name = volume->volumeName();
    if (name.find("MultiLayer") != std::string::npos) {
      nMultiLayers++;
    }
  });

  // we expect to have 24 volumes in the geometry (3 for every sector) and 48
  // multilayers
  BOOST_CHECK_EQUAL(nVolumes, 24u);
  BOOST_CHECK_EQUAL(nMultiLayers, 48u);

  BOOST_REQUIRE(trackingGeometry);
  BOOST_CHECK(trackingGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen3);
  // we expect 5 volumes registered in the tracking geometry highest volume in
  // the hierarchy
  // 3 barrel real volumes and 2 gaps from the cylinder container stack strategy
  BOOST_CHECK_EQUAL(trackingGeometry->highestTrackingVolume()->volumes().size(),
                    5u);

  // visualize the tracking geometry with obj to see that is built correctly
  Acts::ObjVisualization3D obj;
  trackingGeometry->visualize(obj, gctx);
  obj.write("MuonMockupGeometry.obj");
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test
