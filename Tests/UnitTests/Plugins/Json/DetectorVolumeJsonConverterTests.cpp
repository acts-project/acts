// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Geometry/ConeVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "ActsPlugins/Json/DetectorVolumeJsonConverter.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <fstream>
#include <memory>
#include <numbers>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace {
/// Helper method that allows to use the already existing testing
/// infrastructure with the new const-correct detector design
///
std::vector<std::shared_ptr<Surface>> unpackSurfaces(
    const std::vector<Surface*>& surfaces) {
  std::vector<std::shared_ptr<Surface>> uSurfaces;
  uSurfaces.reserve(surfaces.size());
  for (auto* s : surfaces) {
    uSurfaces.push_back(s->getSharedPtr());
  }
  return uSurfaces;
}

}  // namespace

namespace ActsTests {

GeometryContext tContext;
auto cGeometry = CylindricalTrackingGeometry(tContext);

auto portalGenerator = Experimental::defaultPortalGenerator();

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(SingleEmptyVolume) {
  // Create a single cylindrical volume
  Transform3 nominal = Transform3::Identity();
  auto bounds = std::make_unique<CylinderVolumeBounds>(0., 50., 100.);

  auto volume = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "EmptyVolume", nominal, std::move(bounds),
      Experimental::tryAllPortals());

  std::ofstream out;

  auto jVolume =
      DetectorVolumeJsonConverter::toJson(tContext, *volume, {volume.get()});

  out.open("single-empty-volume.json");
  out << jVolume.dump(4);
  out.close();

  auto in = std::ifstream("single-empty-volume.json",
                          std::ifstream::in | std::ifstream::binary);

  BOOST_CHECK(in.good());
  nlohmann::json jVolumeIn;
  in >> jVolumeIn;
  in.close();

  auto volumeIn = DetectorVolumeJsonConverter::fromJson(tContext, jVolumeIn);

  BOOST_CHECK_EQUAL(volumeIn->name(), volume->name());
  BOOST_CHECK(
      volumeIn->transform(tContext).isApprox(volume->transform(tContext)));
  BOOST_CHECK_EQUAL(volumeIn->surfaces().size(), volume->surfaces().size());
  BOOST_CHECK_EQUAL(volumeIn->volumes().size(), volume->volumes().size());
}

BOOST_AUTO_TEST_CASE(SingleSurfaceVolume) {
  // Create a single cylindrical volume
  Transform3 nominal = Transform3::Identity();
  auto volumeBounds = std::make_unique<CylinderVolumeBounds>(0., 50., 100.);
  auto surfaceBounds = std::make_unique<CylinderBounds>(25., 100.);

  auto cylinderSurface =
      Surface::makeShared<CylinderSurface>(nominal, std::move(surfaceBounds));

  auto volume = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "CylinderVolume", nominal,
      std::move(volumeBounds), {cylinderSurface}, {},
      Experimental::tryRootVolumes(), Experimental::tryAllPortalsAndSurfaces());

  std::ofstream out;

  auto jVolume =
      DetectorVolumeJsonConverter::toJson(tContext, *volume, {volume.get()});

  out.open("single-surface-volume.json");
  out << jVolume.dump(4);
  out.close();

  auto in = std::ifstream("single-surface-volume.json",
                          std::ifstream::in | std::ifstream::binary);

  BOOST_CHECK(in.good());
  nlohmann::json jVolumeIn;
  in >> jVolumeIn;
  in.close();

  auto volumeIn = DetectorVolumeJsonConverter::fromJson(tContext, jVolumeIn);

  BOOST_CHECK_EQUAL(volumeIn->name(), volume->name());
  BOOST_CHECK(
      volumeIn->transform(tContext).isApprox(volume->transform(tContext)));
  BOOST_CHECK_EQUAL(volumeIn->surfaces().size(), volume->surfaces().size());
  BOOST_CHECK_EQUAL(volumeIn->volumes().size(), volume->volumes().size());
}

BOOST_AUTO_TEST_CASE(EndcapVolumeWithSurfaces) {
  CylindricalTrackingGeometry::DetectorStore dStore;

  auto rSurfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 36., 0.125, 0.,
                                          55., -800, 2., 22u);

  auto endcapSurfaces =
      std::make_shared<Experimental::LayerStructureBuilder::SurfacesHolder>(
          unpackSurfaces(rSurfaces));
  // Configure the layer structure builder
  Experimental::LayerStructureBuilder::Config lsConfig;
  lsConfig.auxiliary = "*** Endcap with 22 surfaces ***";
  lsConfig.surfacesProvider = endcapSurfaces;
  lsConfig.binnings = {
      {DirectedProtoAxis(AxisDirection::AxisPhi, AxisBoundaryType::Closed,
                         -std::numbers::pi, std::numbers::pi, 22u),
       1u}};

  auto layerBuilder = std::make_shared<Experimental::LayerStructureBuilder>(
      lsConfig, getDefaultLogger("EndcapInteralsBuilder", Logging::VERBOSE));

  Experimental::VolumeStructureBuilder::Config shapeConfig;
  shapeConfig.boundValues = {10, 100, 10., std::numbers::pi, 0.};
  shapeConfig.transform =
      Transform3{Transform3::Identity()}.pretranslate(Vector3(0., 0., -800.));
  shapeConfig.boundsType = VolumeBounds::BoundsType::eCylinder;

  auto shapeBuilder = std::make_shared<Experimental::VolumeStructureBuilder>(
      shapeConfig, getDefaultLogger("EndcapShapeBuilder", Logging::VERBOSE));

  Experimental::DetectorVolumeBuilder::Config dvCfg;
  dvCfg.auxiliary = "*** Test 1 - Cylinder with internal Surface ***";
  dvCfg.name = "CylinderWithSurface";
  dvCfg.externalsBuilder = shapeBuilder;
  dvCfg.internalsBuilder = layerBuilder;

  auto dvBuilder = std::make_shared<Experimental::DetectorVolumeBuilder>(
      dvCfg, getDefaultLogger("EndcapBuilder", Logging::VERBOSE));

  auto [volumes, portals, roots] = dvBuilder->construct(tContext);
  auto volume = volumes.front();

  auto jVolume =
      DetectorVolumeJsonConverter::toJson(tContext, *volume, {volume.get()});

  std::ofstream out;
  out.open("endcap-volume-with-surfaces.json");
  out << jVolume.dump(4);
  out.close();

  auto in = std::ifstream("endcap-volume-with-surfaces.json",
                          std::ifstream::in | std::ifstream::binary);

  BOOST_CHECK(in.good());
  nlohmann::json jVolumeIn;
  in >> jVolumeIn;
  in.close();

  auto volumeIn = DetectorVolumeJsonConverter::fromJson(tContext, jVolumeIn);

  BOOST_CHECK_EQUAL(volumeIn->name(), volume->name());
  BOOST_CHECK(
      volumeIn->transform(tContext).isApprox(volume->transform(tContext)));
  BOOST_CHECK_EQUAL(volumeIn->surfaces().size(), volume->surfaces().size());
  BOOST_CHECK_EQUAL(volumeIn->volumes().size(), volume->volumes().size());

  // Cross-check writing
  jVolume = DetectorVolumeJsonConverter::toJson(tContext, *volumeIn,
                                                {volumeIn.get()});
  out.open("endcap-volume-with-surfaces-closure.json");
  out << jVolume.dump(4);
  out.close();
}

BOOST_AUTO_TEST_CASE(BarrelVolumeWithSurfaces) {
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto cSurfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.145, 72,
                                              3., 2., {32u, 14u});

  auto barrelSurfaces =
      std::make_shared<Experimental::LayerStructureBuilder::SurfacesHolder>(
          unpackSurfaces(cSurfaces));

  // Configure the layer structure builder
  Experimental::LayerStructureBuilder::Config lsConfig;
  lsConfig.auxiliary = "*** Barrel with 448 surfaces ***";
  lsConfig.surfacesProvider = barrelSurfaces;
  lsConfig.binnings = {
      {DirectedProtoAxis(AxisDirection::AxisZ, AxisBoundaryType::Bound, -480.,
                         480., 14u),
       1u},
      {DirectedProtoAxis(AxisDirection::AxisPhi, AxisBoundaryType::Closed,
                         -std::numbers::pi, std::numbers::pi, 32u),
       1u}};

  auto barrelBuilder = std::make_shared<Experimental::LayerStructureBuilder>(
      lsConfig, getDefaultLogger("BarrelInternalsBuilder", Logging::VERBOSE));

  Experimental::VolumeStructureBuilder::Config shapeConfig;
  shapeConfig.boundValues = {60., 80., 800., std::numbers::pi, 0.};
  shapeConfig.boundsType = VolumeBounds::BoundsType::eCylinder;

  auto shapeBuilder = std::make_shared<Experimental::VolumeStructureBuilder>(
      shapeConfig, getDefaultLogger("BarrelShapeBuilder", Logging::VERBOSE));

  Experimental::DetectorVolumeBuilder::Config dvCfg;
  dvCfg.auxiliary = "*** Test 1 - Cylinder with internal Surface ***";
  dvCfg.name = "BarrelWithSurfaces";
  dvCfg.externalsBuilder = shapeBuilder;
  dvCfg.internalsBuilder = barrelBuilder;

  auto dvBuilder = std::make_shared<Experimental::DetectorVolumeBuilder>(
      dvCfg, getDefaultLogger("EndcapBuilder", Logging::VERBOSE));

  auto [volumes, portals, roots] = dvBuilder->construct(tContext);

  auto volume = volumes.front();

  auto jVolume =
      DetectorVolumeJsonConverter::toJson(tContext, *volume, {volume.get()});

  std::ofstream out;
  out.open("barrel-volume-with-surfaces.json");
  out << jVolume.dump(4);
  out.close();

  auto in = std::ifstream("barrel-volume-with-surfaces.json",
                          std::ifstream::in | std::ifstream::binary);

  BOOST_CHECK(in.good());
  nlohmann::json jVolumeIn;
  in >> jVolumeIn;
  in.close();

  auto volumeIn = DetectorVolumeJsonConverter::fromJson(tContext, jVolumeIn);

  BOOST_CHECK_EQUAL(volumeIn->name(), volume->name());
  BOOST_CHECK(
      volumeIn->transform(tContext).isApprox(volume->transform(tContext)));
  BOOST_CHECK_EQUAL(volumeIn->surfaces().size(), volume->surfaces().size());
  BOOST_CHECK_EQUAL(volumeIn->volumes().size(), volume->volumes().size());

  // Cross-check writing
  jVolume = DetectorVolumeJsonConverter::toJson(tContext, *volumeIn,
                                                {volumeIn.get()});
  out.open("barrel-volume-with-surfaces-closure.json");
  out << jVolume.dump(4);
  out.close();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
