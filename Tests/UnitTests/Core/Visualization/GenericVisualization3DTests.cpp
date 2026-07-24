// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Visualization/GenericVisualization3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <numbers>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace Acts;

namespace ActsTests {

namespace {

/// The relevant records of an obj stream: vertex positions, 0-based face
/// index lists and 0-based line index pairs
struct ParsedObj {
  std::vector<Vector3> vertices;
  std::vector<std::vector<std::size_t>> faces;
  std::vector<std::pair<std::size_t, std::size_t>> lines;
};

/// Parse the `v`, `f` and `l` records from an obj stream
/// @param objString the obj stream content
/// @return the parsed records
ParsedObj parseObj(const std::string& objString) {
  ParsedObj parsed;
  auto ss = std::stringstream{objString};
  for (std::string line; std::getline(ss, line);) {
    auto ls = std::stringstream{line};
    std::string tag;
    ls >> tag;
    if (tag == "v") {
      Vector3 vtx = Vector3::Zero();
      ls >> vtx.x() >> vtx.y() >> vtx.z();
      parsed.vertices.push_back(vtx);
    } else if (tag == "f") {
      std::vector<std::size_t> face;
      for (std::size_t idx = 0; ls >> idx;) {
        face.push_back(idx - 1);
      }
      parsed.faces.push_back(std::move(face));
    } else if (tag == "l") {
      std::size_t start = 0;
      std::size_t end = 0;
      ls >> start >> end;
      parsed.lines.push_back({start - 1, end - 1});
    }
  }
  return parsed;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(VisualizationSuite)

BOOST_AUTO_TEST_CASE(GenericVisualization3DVertexLine) {
  GenericVisualization3D helper;

  helper.vertex({1., 2., 3.}, {10, 20, 30});
  BOOST_REQUIRE_EQUAL(helper.vertices().size(), 1);
  CHECK_CLOSE_ABS(helper.vertices()[0].position, Vector3(1., 2., 3.), 1e-9);
  BOOST_CHECK_EQUAL(helper.vertices()[0].color, Color(10, 20, 30));
  BOOST_CHECK(helper.faces().empty());
  BOOST_CHECK(helper.lines().empty());

  helper.line({0., 0., 1.}, {1., 0., 0.}, {40, 50, 60});
  BOOST_REQUIRE_EQUAL(helper.lines().size(), 1);
  CHECK_CLOSE_ABS(helper.lines()[0].a, Vector3(0., 0., 1.), 1e-9);
  CHECK_CLOSE_ABS(helper.lines()[0].b, Vector3(1., 0., 0.), 1e-9);
  BOOST_CHECK_EQUAL(helper.lines()[0].color, Color(40, 50, 60));
  // The line end points are not appended to the vertex collection
  BOOST_CHECK_EQUAL(helper.vertices().size(), 1);

  helper.clear();
  BOOST_CHECK(helper.vertices().empty());
  BOOST_CHECK(helper.faces().empty());
  BOOST_CHECK(helper.lines().empty());
}

BOOST_AUTO_TEST_CASE(GenericVisualization3DFaces) {
  GenericVisualization3D helper;

  // face() appends the corner points as vertices and one indexed face
  const std::vector<Vector3> triangle = {
      {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.}};
  helper.face(triangle, {10, 20, 30});
  BOOST_REQUIRE_EQUAL(helper.vertices().size(), 3);
  BOOST_REQUIRE_EQUAL(helper.faces().size(), 1);
  BOOST_CHECK(helper.faces()[0].indices ==
              (IVisualization3D::FaceType{0, 1, 2}));
  BOOST_CHECK_EQUAL(helper.faces()[0].color, Color(10, 20, 30));
  for (std::size_t i = 0; i < triangle.size(); ++i) {
    CHECK_CLOSE_ABS(helper.vertices()[i].position, triangle[i], 1e-9);
    BOOST_CHECK_EQUAL(helper.vertices()[i].color, Color(10, 20, 30));
  }

  // faces() with an empty face list behaves like face()
  GenericVisualization3D viaFaces;
  viaFaces.faces(triangle, {}, {10, 20, 30});
  BOOST_CHECK_EQUAL(viaFaces.vertices().size(), helper.vertices().size());
  BOOST_REQUIRE_EQUAL(viaFaces.faces().size(), helper.faces().size());
  BOOST_CHECK(viaFaces.faces()[0].indices == helper.faces()[0].indices);

  // faces() offsets the indices by the previously collected vertices,
  // two-index faces are stored as lines from the referenced points
  const std::vector<Vector3> quad = {
      {0., 0., 0.}, {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.}};
  helper.faces(quad, {{0, 1, 2, 3}, {0, 2}}, {40, 50, 60});
  BOOST_REQUIRE_EQUAL(helper.vertices().size(), 7);
  BOOST_REQUIRE_EQUAL(helper.faces().size(), 2);
  BOOST_CHECK(helper.faces()[1].indices ==
              (IVisualization3D::FaceType{3, 4, 5, 6}));
  BOOST_CHECK_EQUAL(helper.faces()[1].color, Color(40, 50, 60));
  BOOST_REQUIRE_EQUAL(helper.lines().size(), 1);
  CHECK_CLOSE_ABS(helper.lines()[0].a, quad[0], 1e-9);
  CHECK_CLOSE_ABS(helper.lines()[0].b, quad[2], 1e-9);
}

BOOST_AUTO_TEST_CASE(GenericVisualization3DConstruction) {
  // Conformance to the IVisualization3D interface, cf.
  // Visualization3DConstruction in Visualization3DTests.cpp
  GenericVisualization3D generic;
  IVisualization3D* vis = &generic;
  vis->vertex({0., 0., 0.});
  vis->line({0., 0., 0.}, {1., 0., 0.});

  std::stringstream ss;
  ss << *vis;
  BOOST_CHECK_EQUAL(ss.str(),
                    "GenericVisualization3D: 1 vertices, 0 faces, 1 lines\n");

  vis->clear();
  BOOST_CHECK(generic.vertices().empty());
  BOOST_CHECK(generic.lines().empty());
}

/// Draw the same surfaces with ObjVisualization3D and GenericVisualization3D
/// and check that the collected data matches the obj output element-for-element
BOOST_AUTO_TEST_CASE(GenericVisualization3DAgainstObj) {
  auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto identity = Transform3::Identity();

  const double halfPhiSector = std::numbers::pi / 4.;

  // The reference surfaces, same parameters as in SurfaceView3DBase.hpp
  std::vector<std::shared_ptr<Surface>> surfaces;
  surfaces.push_back(Surface::makeShared<PlaneSurface>(
      identity, std::make_shared<RectangleBounds>(2., 3.)));
  surfaces.push_back(Surface::makeShared<CylinderSurface>(
      identity, std::make_shared<CylinderBounds>(5., 10.)));
  surfaces.push_back(Surface::makeShared<ConeSurface>(
      identity, std::make_shared<ConeBounds>(0.245, 0., 10., halfPhiSector)));

  for (const auto& surface : surfaces) {
    // High output precision for the text based comparison
    ObjVisualization3D obj{10};
    GenericVisualization3D generic;
    GeometryView3D::drawSurface(obj, *surface, gctx, identity, s_viewSensitive);
    GeometryView3D::drawSurface(generic, *surface, gctx, identity,
                                s_viewSensitive);

    std::stringstream objStream;
    obj.write(objStream);
    const ParsedObj parsed = parseObj(objStream.str());

    BOOST_REQUIRE(!generic.vertices().empty());
    BOOST_REQUIRE_EQUAL(generic.vertices().size(), parsed.vertices.size());
    for (std::size_t i = 0; i < parsed.vertices.size(); ++i) {
      CHECK_CLOSE_ABS(generic.vertices()[i].position, parsed.vertices[i], 1e-6);
      BOOST_CHECK_EQUAL(generic.vertices()[i].color, s_viewSensitive.color);
    }

    BOOST_REQUIRE(!generic.faces().empty());
    BOOST_REQUIRE_EQUAL(generic.faces().size(), parsed.faces.size());
    for (std::size_t i = 0; i < parsed.faces.size(); ++i) {
      BOOST_CHECK(generic.faces()[i].indices == parsed.faces[i]);
      BOOST_CHECK_EQUAL(generic.faces()[i].color, s_viewSensitive.color);
    }

    // Lines reference obj vertices, compare the end points
    BOOST_REQUIRE_EQUAL(generic.lines().size(), parsed.lines.size());
    for (std::size_t i = 0; i < parsed.lines.size(); ++i) {
      CHECK_CLOSE_ABS(generic.lines()[i].a,
                      parsed.vertices[parsed.lines[i].first], 1e-6);
      CHECK_CLOSE_ABS(generic.lines()[i].b,
                      parsed.vertices[parsed.lines[i].second], 1e-6);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
