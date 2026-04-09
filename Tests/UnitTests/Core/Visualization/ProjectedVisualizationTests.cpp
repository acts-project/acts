// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ProjectedVisualization.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(VisualizationSuite)

BOOST_AUTO_TEST_CASE(ProjectedVisualizationConstruction) {
  ProjectedVisualization vis({
      {"xy", projectToXY},
      {"zphi", projectToZPhi},
      {"zr", projectToZR},
  });

  BOOST_CHECK_EQUAL(vis.projections().size(), 3u);
  BOOST_CHECK(vis.projections().contains("xy"));
  BOOST_CHECK(vis.projections().contains("zphi"));
  BOOST_CHECK(vis.projections().contains("zr"));

  IVisualization3D* base = &vis;
  std::stringstream ss;
  base->write(ss);

  const std::string out = ss.str();
  BOOST_CHECK(out.find("Projection: xy") != std::string::npos);
  BOOST_CHECK(out.find("Projection: zphi") != std::string::npos);
  BOOST_CHECK(out.find("Projection: zr") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(ProjectedVisualizationProjectionFunctions) {
  const Vector3 p{3., 4., 5.};

  const auto xy = projectToXY(p);
  CHECK_CLOSE_ABS(xy.x(), 3., 1e-12);
  CHECK_CLOSE_ABS(xy.y(), 4., 1e-12);

  const auto zphi = projectToZPhi(p);
  CHECK_CLOSE_ABS(zphi.x(), 5., 1e-12);
  CHECK_CLOSE_ABS(zphi.y(), std::atan2(4., 3.), 1e-12);

  const auto zr = projectToZR(p);
  CHECK_CLOSE_ABS(zr.x(), 5., 1e-12);
  CHECK_CLOSE_ABS(zr.y(), 5., 1e-12);
}

BOOST_AUTO_TEST_CASE(ProjectedVisualizationWriteAndClear) {
  ProjectedVisualization vis({{"xy", projectToXY}});

  vis.vertex({1., 2., 3.}, {1, 2, 3});
  vis.line({0., 0., 0.}, {1., 1., 1.}, {4, 5, 6});
  vis.face({{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}}, {7, 8, 9});
  vis.faces({{0., 0., 0.}, {2., 0., 0.}, {2., 2., 0.}, {0., 2., 0.}},
            {{0, 1, 2}, {0, 2, 3}}, {11, 12, 13});

  const auto& projectedVertices = vis.projectedVertices();
  BOOST_CHECK(projectedVertices.contains("xy"));
  BOOST_CHECK_EQUAL(projectedVertices.at("xy").size(), 1u);
  const auto& [vtx, vtxColor] = projectedVertices.at("xy").front();
  CHECK_CLOSE_ABS(vtx.x(), 1., 1e-12);
  CHECK_CLOSE_ABS(vtx.y(), 2., 1e-12);
  BOOST_CHECK_EQUAL(vtxColor[0], 1);
  BOOST_CHECK_EQUAL(vtxColor[1], 2);
  BOOST_CHECK_EQUAL(vtxColor[2], 3);

  const auto& projectedLines = vis.projectedLines();
  BOOST_CHECK(projectedLines.contains("xy"));
  BOOST_CHECK_EQUAL(projectedLines.at("xy").size(), 1u);
  const auto& [line, lineColor] = projectedLines.at("xy").front();
  CHECK_CLOSE_ABS(line[0].x(), 0., 1e-12);
  CHECK_CLOSE_ABS(line[0].y(), 0., 1e-12);
  CHECK_CLOSE_ABS(line[1].x(), 1., 1e-12);
  CHECK_CLOSE_ABS(line[1].y(), 1., 1e-12);
  BOOST_CHECK_EQUAL(lineColor[0], 4);
  BOOST_CHECK_EQUAL(lineColor[1], 5);
  BOOST_CHECK_EQUAL(lineColor[2], 6);

  const auto& projectedFaces = vis.projectedFaces();
  BOOST_CHECK(projectedFaces.contains("xy"));
  BOOST_CHECK_EQUAL(projectedFaces.at("xy").size(), 3u);

  std::stringstream ss;
  vis.write(ss);
  const std::string out = ss.str();

  BOOST_CHECK(out.find("Projection: xy") != std::string::npos);
  BOOST_CHECK(out.find("Vertex: 1 2 1 2 3") != std::string::npos);
  BOOST_CHECK(out.find("Line: 0 0 1 1 4 5 6") != std::string::npos);
  BOOST_CHECK(out.find("Face: 3 7 8 9") != std::string::npos);

  vis.clear();
  BOOST_CHECK(vis.projectedVertices().empty());
  BOOST_CHECK(vis.projectedLines().empty());
  BOOST_CHECK(vis.projectedFaces().empty());

  std::stringstream cleared;
  vis.write(cleared);
  const std::string clearedOut = cleared.str();

  BOOST_CHECK(clearedOut.find("Projection: xy") != std::string::npos);
  BOOST_CHECK(clearedOut.find("Vertex:") == std::string::npos);
  BOOST_CHECK(clearedOut.find("Line:") == std::string::npos);
  BOOST_CHECK(clearedOut.find("Face:") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(ProjectedVisualizationWriteToPath) {
  ProjectedVisualization vis({{"xy", projectToXY}});
  vis.vertex({9., 8., 7.}, {10, 11, 12});

  const auto path =
      std::filesystem::temp_directory_path() /
      std::filesystem::path("acts_projected_visualization_test.txt");

  vis.write(path);

  std::ifstream ifs(path.string());
  BOOST_CHECK(ifs.good());

  const std::string content((std::istreambuf_iterator<char>(ifs)),
                            std::istreambuf_iterator<char>());
  BOOST_CHECK(content.find("Projection: xy") != std::string::npos);
  BOOST_CHECK(content.find("Vertex: 9 8 10 11 12") != std::string::npos);

  std::filesystem::remove(path);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
