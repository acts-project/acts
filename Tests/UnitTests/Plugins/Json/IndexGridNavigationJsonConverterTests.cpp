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
#include "ActsPlugins/Json/IndexGridNavigationJsonConverter.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <nlohmann/json.hpp>
#include <fstream>

using namespace Acts;

namespace ActsTests {

GeometryContext tContext;
CylindricalTrackingGeometry cGeometry = CylindricalTrackingGeometry(tContext);
CylindricalTrackingGeometry::DetectorStore dStore;

auto tLogger =
    getDefaultLogger("IndexGridNavigationJsonConverterTests", Logging::VERBOSE);

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(RegularCylinderIndexGridNavigationPolicyToJson) {
  ACTS_LOCAL_LOGGER(getDefaultLogger(
      "*** Test IndexGridNavigationJsonConverter", Logging::VERBOSE));
  ACTS_INFO(
      "Testing RegularCylinderIndexGridNavigationPolicy to Json conversion.");

  // Get surfaces and dump them into a Tracking volume
  auto bSurfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.14, 31.,
                                              3., 2., {16, 14});

  TrackingVolume tVolume(
      Transform3::Identity(),
      std::make_shared<Acts::CylinderVolumeBounds>(0., 50., 500.),
      "CylinderVolume");


  nlohmann::json jSurfaceVertices = nlohmann::json::array();    
  // Measure the grid extents
  Acts::Extent extent;
  for (const auto& [idx, surface] : Acts::enumerate(bSurfaces)) {
    auto phSurface = surface->polyhedronRepresentation(tContext, 1);
    std::vector<std::array<double,3>> vertex;
    for (const auto& vtx : phSurface.vertices){
        vertex.push_back({vtx.x(), vtx.y(), vtx.z()});
    }
    jSurfaceVertices.push_back(vertex);
    extent.extend(phSurface.extent());
    tVolume.addSurface(surface->getSharedPtr());

  }

  // z-phi Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound>  axisZ(extent.min(AxisDirection::AxisZ), extent.max(AxisDirection::AxisZ),
             14);
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(
      -std::numbers::pi, std::numbers::pi, 16);
  Grid gridZPhi(Type<std::vector<std::size_t>>, std::move(axisZ),
                std::move(axisPhi));

  // Indexed Surface grid
  IndexGrid<decltype(gridZPhi)> indexedGridZPhi(
      std::move(gridZPhi), {AxisDirection::AxisZ, AxisDirection::AxisPhi});

  Acts::Experimental::IndexGridNavigationConfig igCfg;
  igCfg.binExpansion = {0u, 0u};
  //igCfg.referenceGenerator = std::make_shared<CenterReferenceGenerator>();
  igCfg.referenceGenerator = std::make_shared<PolyhedronReferenceGenerator>();

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

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests