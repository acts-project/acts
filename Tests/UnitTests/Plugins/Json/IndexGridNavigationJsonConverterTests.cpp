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
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/IndexGridNavigationPolicy.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsPlugins/Json/IndexGridNavigationJsonConverter.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <fstream>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

GeometryContext tContext;
CylindricalTrackingGeometry cGeometry = CylindricalTrackingGeometry(tContext);
CylindricalTrackingGeometry::DetectorStore dStore;

auto tLogger =
    getDefaultLogger("IndexGridNavigationJsonConverterTests", Logging::VERBOSE);

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(RegularCylinderIndexGridNavigationWithProjection) {
  double ringRadius = 36.;

  auto referenceSurface = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(),
      std::make_shared<CylinderBounds>(ringRadius, 10e10));

  // The Projected reference geneator wihtout expansion
  auto projRfGen = std::make_shared<ProjectedReferenceGenerator>();
  projRfGen->referenceSurface = referenceSurface;
  projRfGen->luminousRegion = { Vector3{0.,0.,0.} }; //{Vector3{0., 0., -150.}, Vector3{0., 0., 150.}};

  double phi = 0.;
  double ringZ = 0.;
  double surfaceTilt = 0.0;

  Vector3 surfaceCenter(ringRadius * std::cos(phi), ringRadius * std::sin(phi),
                        ringZ);

  // Local z axis is the normal vector
  Vector3 surfaceLocalZ(std::cos(phi + surfaceTilt),
                        std::sin(phi + surfaceTilt), 0.);
  // Local y axis is the global z axis
  Vector3 surfaceLocalY(0., 0., 1);
  // Local x axis the normal to local y,z
  Vector3 surfaceLocalX(-std::sin(phi + surfaceTilt),
                        std::cos(phi + surfaceTilt), 0.);
  // Create the RotationMatrix
  RotationMatrix3 surfaceRotation;
  surfaceRotation.col(0) = surfaceLocalX;
  surfaceRotation.col(1) = surfaceLocalY;
  surfaceRotation.col(2) = surfaceLocalZ;
  // Get the surfaceTransform
  auto surfaceTransform =
      Transform3(Translation3(surfaceCenter) * surfaceRotation);

  auto centralSurface = Surface::makeShared<PlaneSurface>(
      surfaceTransform, std::make_shared<RectangleBounds>(12., 20.));

  // Grap a central surface in z and phi
  // auto centralSurface =

  nlohmann::json jSurfaceVertices = nlohmann::json::array();
  auto phSurface = centralSurface->polyhedronRepresentation(tContext, 1);
  std::vector<std::array<double, 3>> vertex;
  for (const auto& vtx : phSurface.vertices) {
    vertex.push_back({vtx.x(), vtx.y(), vtx.z()});
  }
  jSurfaceVertices.push_back(vertex);

  // z-phi Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisZ(-100, 100, 50);
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(-1., 1., 10);
  Grid gridZPhi(Type<std::vector<std::size_t>>, std::move(axisZ),
                std::move(axisPhi));

  TrackingVolume tVolume(Transform3::Identity(),
                         std::make_shared<CylinderVolumeBounds>(0., 50., 50.),
                         "CylinderVolume");

  tVolume.addSurface(centralSurface->getSharedPtr());

  // Indexed Surface grid
  IndexGrid<decltype(gridZPhi)> indexedGridZPhi(
      std::move(gridZPhi), {AxisDirection::AxisZ, AxisDirection::AxisPhi});

  Acts::Experimental::IndexGridNavigationConfig igCfg;
  igCfg.binExpansion = {0u, 0u};
  igCfg.referenceExpansion = {0., 0.15};
  igCfg.referenceGenerator = projRfGen;
  Experimental::RegularCylinderIndexGridNavigationPolicy policy(
      tContext, tVolume, *tLogger, igCfg, indexedGridZPhi);

  nlohmann::json jPolicy =
      Acts::IndexGridNavigationJsonConverter::toJson(policy);

  // Output surface vertices
  nlohmann::json outputJson;

  outputJson["SurfaceVertices"] = jSurfaceVertices;
  outputJson["NavigationPolicy"] = jPolicy;

  // Write to file
  std::ofstream jsonFile("RegularCylinderIndexGridNavigationPolicy.json");
  jsonFile << std::setw(2) << outputJson;
  jsonFile.close();
}

BOOST_AUTO_TEST_CASE(RegularCylinderIndexGridNavigationPolicyToJson) {
  ACTS_LOCAL_LOGGER(getDefaultLogger(
      "*** Test IndexGridNavigationJsonConverter", Logging::VERBOSE));
  ACTS_INFO(
      "Testing RegularCylinderIndexGridNavigationPolicy to Json conversion.");

  // Get surfaces and dump them into a Tracking volume
  auto bSurfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.14, 31.,
                                              3., 2., {16, 14});

  // Pre-loop for getting the extent
  Extent extent;
  for (const auto surface : bSurfaces) {
    auto phSurface = surface->polyhedronRepresentation(tContext, 1);
    extent.extend(phSurface.extent());
  }
  auto referenceSurface = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(
                                  extent.medium(AxisDirection::AxisR), 10e10));

  // The Projected reference geneator wihtout expansion
  auto projRfGen = std::make_shared<ProjectedReferenceGenerator>();
  projRfGen->referenceSurface = referenceSurface;
  projRfGen->luminousRegion = {Vector3{0., 0., -150.}, Vector3{0., 0., 150.}};

  // We test this for different Reference generators
  std::vector<std::tuple<std::string, std::array<std::size_t, 2>,
                         std::shared_ptr<IReferenceGenerator>>>
      referenceGenerators = {
          {"Polyhedron_b_14_16",
           {14, 16},
           std::make_shared<PolyhedronReferenceGenerator>()},
          {"Polyhedron_b_14_16_m3",
           {14 * 3, 16 * 3},
           std::make_shared<PolyhedronReferenceGenerator>()},
          {"Projected_b_14_16", {14, 16}, projRfGen},
          {"Projected_b_14_16_m3", {14 * 3, 16 * 3}, projRfGen},
          {"CenterReferenceGenerator",
           {14, 16},
           std::make_shared<CenterReferenceGenerator>()}};

  for (auto [name, bins, refGenerator] : referenceGenerators) {
    ACTS_INFO(" Testing Reference Generator: " << name);

    TrackingVolume tVolume(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 50., 500.),
        "CylinderVolume");

    nlohmann::json jSurfaceVertices = nlohmann::json::array();
    // Measure the grid extents
    for (const auto& [idx, surface] : enumerate(bSurfaces)) {
      auto phSurface = surface->polyhedronRepresentation(tContext, 1);
      std::vector<std::array<double, 3>> vertex;
      for (const auto& vtx : phSurface.vertices) {
        vertex.push_back({vtx.x(), vtx.y(), vtx.z()});
      }
      jSurfaceVertices.push_back(vertex);
      tVolume.addSurface(surface->getSharedPtr());
    }

    // z-phi Axes & Grid
    Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisZ(
        extent.min(AxisDirection::AxisZ), extent.max(AxisDirection::AxisZ),
        bins[0]);
    Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(
        -std::numbers::pi, std::numbers::pi, bins[1]);
    Grid gridZPhi(Type<std::vector<std::size_t>>, std::move(axisZ),
                  std::move(axisPhi));

    // Indexed Surface grid
    IndexGrid<decltype(gridZPhi)> indexedGridZPhi(
        std::move(gridZPhi), {AxisDirection::AxisZ, AxisDirection::AxisPhi});

    Acts::Experimental::IndexGridNavigationConfig igCfg;
    igCfg.binExpansion = {0u, 0u};
    if (name != "CenterReferenceGenerator") {
      igCfg.referenceExpansion = {10., 0.0};
    }

    // igCfg.referenceGenerator = std::make_shared<CenterReferenceGenerator>();
    igCfg.referenceGenerator = refGenerator;
    Experimental::RegularCylinderIndexGridNavigationPolicy policy(
        tContext, tVolume, *tLogger, igCfg, indexedGridZPhi);

    nlohmann::json jPolicy =
        Acts::IndexGridNavigationJsonConverter::toJson(policy);

    // Output surface vertices
    nlohmann::json outputJson;

    outputJson["SurfaceVertices"] = jSurfaceVertices;
    outputJson["NavigationPolicy"] = jPolicy;

    // Write to file
    std::ofstream jsonFile("RegularCylinderIndexGridNavigationPolicy_" + name +
                           ".json");
    jsonFile << std::setw(2) << outputJson;
    jsonFile.close();
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests