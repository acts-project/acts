// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/IndexGridNavigationPolicy.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "ActsPlugins/Json/IndexGridNavigationJsonConverter.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <fstream>
#include <numbers>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

auto tContext = GeometryContext::dangerouslyDefaultConstruct();
CylindricalTrackingGeometry cGeometry = CylindricalTrackingGeometry(tContext);
CylindricalTrackingGeometry::DetectorStore dStore;

auto tLogger =
    getDefaultLogger("IndexGridNavigationJsonConverterTests", Logging::INFO);

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(RegularCylinderIndexGridNavigationPolicyToJson) {
  auto tContext = GeometryContext::dangerouslyDefaultConstruct();
  CylindricalTrackingGeometry cGeometry = CylindricalTrackingGeometry(tContext);
  CylindricalTrackingGeometry::DetectorStore dStore;

  // Get surfaces and dump them into a Tracking volume
  auto bSurfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.14, 31.,
                                              3., 2., {16, 14});

  TrackingVolume tVolume(Transform3::Identity(),
                         std::make_shared<CylinderVolumeBounds>(0., 50., 500.),
                         "CylinderVolume");

  for (const auto surface : bSurfaces) {
    tVolume.addSurface(surface->getSharedPtr());
  }

  // Create the phi-z axes and grids
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(
      -std::numbers::pi, std::numbers::pi, 36);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisZ(-400., 400., 30);
  Grid gridPhiZ(Type<std::vector<std::size_t>>, std::move(axisPhi),
                std::move(axisZ));

  IndexGrid<decltype(gridPhiZ)> indexedGridPhiZ(
      std::move(gridPhiZ), {AxisDirection::AxisPhi, AxisDirection::AxisZ});
  IndexGridNavigationConfig phConfig;
  phConfig.referenceGenerator =
      std::make_shared<PolyhedronReferenceGenerator>();
  IndexGridNavigationPolicy<decltype(gridPhiZ)> polyNavigationPolicy(
      tContext, tVolume, *tLogger, phConfig, indexedGridPhiZ);

  // Convert to Json
  IndexGridNavigationJsonConverter::Options options;
  options.writeSurfaces = true;
  options.writeProjectedSurfaces = true;
  nlohmann::json outputJson = IndexGridNavigationJsonConverter::toJson(
      tContext, polyNavigationPolicy, tVolume, options);
  BOOST_CHECK(outputJson.contains("surfaces"));
  BOOST_CHECK(outputJson.contains("projectedSurfaces"));
  BOOST_CHECK(outputJson.contains("type"));
  BOOST_CHECK(outputJson.contains("grid"));
  BOOST_CHECK_EQUAL(outputJson["type"],
                    "RegularCylinderIndexGridNavigationPolicy");
  // Write to file
  std::ofstream jsonFile("RegularCylinderIndexGridNavigationPolicy.json");
  jsonFile << std::setw(2) << outputJson;
  jsonFile.close();
}

BOOST_AUTO_TEST_CASE(RegularDiscIndexGridNavigationPolicyToJson) {
  TrackingVolume tVolume(Transform3::Identity(),
                         std::make_shared<CylinderVolumeBounds>(0., 10., 2.),
                         "DiscVolume");

  // Let's make a few disc segment surfaces
  std::vector<std::shared_ptr<Surface>> discSurfaces;
  double rIr = 3.;
  auto boundsIr = std::make_shared<TrapezoidBounds>(0.2, 1.4, 2.2);

  double rOr = 6.;
  auto boundsOr = std::make_shared<TrapezoidBounds>(1.2, 2.3, 2.2);

  for (std::size_t is = 0; is < 12; ++is) {
    double phiIr = -std::numbers::pi + is * (std::numbers::pi / 6.0);
    AngleAxis3 rotationIr(phiIr - 0.5 * std::numbers::pi, Vector3::UnitZ());
    Transform3 transformIr(rotationIr);
    transformIr.translation() =
        Vector3(rIr * std::cos(phiIr), rIr * std::sin(phiIr), 0.);
    tVolume.addSurface(
        Surface::makeShared<PlaneSurface>(transformIr, boundsIr));

    double phiOr = -std::numbers::pi + (is + 1) * (std::numbers::pi / 6.0);
    AngleAxis3 rotationOr(phiOr - 0.5 * std::numbers::pi, Vector3::UnitZ());
    Transform3 transformOr(rotationOr);
    transformOr.translation() =
        Vector3(rOr * std::cos(phiOr), rOr * std::sin(phiOr), 0.);
    tVolume.addSurface(
        Surface::makeShared<PlaneSurface>(transformOr, boundsOr));
  }

  // Create the r-phi axes and grids
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisR(0.6, 8.6, 8);
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(
      -std::numbers::pi, std::numbers::pi, 36);
  Grid gridRPhi(Type<std::vector<std::size_t>>, std::move(axisR),
                std::move(axisPhi));

  IndexGrid<decltype(gridRPhi)> indexedGridRPhi(
      std::move(gridRPhi), {AxisDirection::AxisR, AxisDirection::AxisPhi});

  IndexGridNavigationConfig phConfig;
  auto referenceGenerator = std::make_shared<PolyhedronReferenceGenerator>();
  referenceGenerator->nSegements = 1;
  phConfig.referenceGenerator = referenceGenerator;

  IndexGridNavigationPolicy<decltype(gridRPhi)> polyNavigationPolicy(
      tContext, tVolume, *tLogger, phConfig, indexedGridRPhi);

  // Convert to Json
  IndexGridNavigationJsonConverter::Options options;
  options.writeSurfaces = true;
  options.writeProjectedSurfaces = true;
  nlohmann::json outputJson = IndexGridNavigationJsonConverter::toJson(
      tContext, polyNavigationPolicy, tVolume, options);

  BOOST_CHECK(outputJson.contains("surfaces"));
  BOOST_CHECK(outputJson.contains("projectedSurfaces"));
  BOOST_CHECK(outputJson.contains("type"));
  BOOST_CHECK(outputJson.contains("grid"));
  BOOST_CHECK_EQUAL(outputJson["type"], "RegularDiscIndexGridNavigationPolicy");

  // Write to file
  std::ofstream jsonFile("RegularDiscIndexGridNavigationPolicy.json");
  jsonFile << std::setw(2) << outputJson;
  jsonFile.close();
}

BOOST_AUTO_TEST_CASE(RegularRingIndexGridNavigationPolicyToJson) {
  // Create a tracking volume
  TrackingVolume tVolume(Transform3::Identity(),
                         std::make_shared<CylinderVolumeBounds>(0., 50., 5.),
                         "RingVolume");
  // Fill with a few disc surfaces
  for (std::size_t is = 0; is < 16; ++is) {
    double phi = -std::numbers::pi + is * (std::numbers::pi / 8.0);
    AngleAxis3 rotation(phi - 0.5 * std::numbers::pi, Vector3::UnitZ());
    Transform3 transform(rotation);
    transform.translation() =
        Vector3(10 * std::cos(phi), 10 * std::sin(phi), 0.);
    tVolume.addSurface(Surface::makeShared<PlaneSurface>(
        transform, std::make_shared<TrapezoidBounds>(1.8, 3.6, 8.0)));
  }

  // Create the r-z axes and grids
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(
      -std::numbers::pi, std::numbers::pi, 8);
  Grid gridPhi(Type<std::vector<std::size_t>>, std::move(axisPhi));

  IndexGrid<decltype(gridPhi)> indexedGridPhi(std::move(gridPhi),
                                              {AxisDirection::AxisPhi});
  IndexGridNavigationConfig phConfig;
  phConfig.referenceGenerator =
      std::make_shared<PolyhedronReferenceGenerator>();
  IndexGridNavigationPolicy<decltype(gridPhi)> polyNavigationPolicy(
      tContext, tVolume, *tLogger, phConfig, indexedGridPhi);

  // Convert to Json
  IndexGridNavigationJsonConverter::Options options;
  options.numPolyhedronPoints = 2;
  options.writeSurfaces = true;
  options.writeProjectedSurfaces = true;
  nlohmann::json outputJson = IndexGridNavigationJsonConverter::toJson(
      tContext, polyNavigationPolicy, tVolume, options);

  BOOST_CHECK(outputJson.contains("surfaces"));
  BOOST_CHECK(outputJson.contains("projectedSurfaces"));
  BOOST_CHECK(outputJson.contains("type"));
  BOOST_CHECK(outputJson.contains("grid"));
  BOOST_CHECK(outputJson.contains("projectedReferenceRange"));
  BOOST_CHECK_EQUAL(outputJson["type"], "RegularRingIndexGridNavigationPolicy");

  // Write to file
  std::ofstream jsonFile("RegularRingIndexGridNavigationPolicy.json");
  jsonFile << std::setw(2) << outputJson;
  jsonFile.close();
}

BOOST_AUTO_TEST_CASE(RegularPlaneIndexGridNavigationPolicyToJson) {
  // Let's create a simple plane in the XY plane
  auto transform0 = Transform3::Identity();
  transform0.translation() = Vector3(-2., -2., 0.);
  auto planeSurface0 = Surface::makeShared<PlaneSurface>(
      transform0, std::make_shared<RectangleBounds>(2.49, 1.99));

  auto transform1 = Transform3::Identity();
  transform1.translation() = Vector3(-1.5, 2., 0.);
  auto planeSurface1 = Surface::makeShared<PlaneSurface>(
      transform1, std::make_shared<RectangleBounds>(2.99, 1.99));

  auto transform2 = Transform3::Identity();
  transform2.translation() = Vector3(2., 0., 0.);
  auto planeSurface2 = Surface::makeShared<PlaneSurface>(
      transform2, std::make_shared<RectangleBounds>(2.49, 3.49));

  // x-y Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(-4.5, 4.5, 9);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(-3.5, 3.5, 7);
  Grid gridXY(Type<std::vector<std::size_t>>, std::move(axisX),
              std::move(axisY));

  TrackingVolume tVolume(Transform3::Identity(),
                         std::make_shared<CuboidVolumeBounds>(20., 20., 5.),
                         "CuboidVolume");
  tVolume.addSurface(planeSurface0);
  tVolume.addSurface(planeSurface1);
  tVolume.addSurface(planeSurface2);

  IndexGrid<decltype(gridXY)> indexedGridXY(
      std::move(gridXY), {AxisDirection::AxisX, AxisDirection::AxisY});
  IndexGridNavigationConfig phConfig;
  phConfig.referenceGenerator =
      std::make_shared<PolyhedronReferenceGenerator>();
  IndexGridNavigationPolicy<decltype(gridXY)> polyNavigationPolicy(
      tContext, tVolume, *tLogger, phConfig, indexedGridXY);
  // Convert to Json
  IndexGridNavigationJsonConverter::Options options;
  options.writeSurfaces = true;
  options.writeProjectedSurfaces = true;
  nlohmann::json outputJson = IndexGridNavigationJsonConverter::toJson(
      tContext, polyNavigationPolicy, tVolume, options);

  BOOST_CHECK(outputJson.contains("surfaces"));
  BOOST_CHECK(outputJson.contains("projectedSurfaces"));
  BOOST_CHECK(outputJson.contains("type"));
  BOOST_CHECK(outputJson.contains("grid"));
  BOOST_CHECK_EQUAL(outputJson["type"],
                    "RegularPlaneIndexGridNavigationPolicy");

  // Write to file
  std::ofstream jsonFile("RegularPlaneIndexGridNavigationPolicy.json");
  jsonFile << std::setw(2) << outputJson;
  jsonFile.close();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
