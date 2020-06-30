// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

// Helper
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <fstream>
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Visualization/GeometryView.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"

namespace Acts {

using namespace UnitLiterals;

namespace Test {

BOOST_AUTO_TEST_SUITE(Geometry)

/// Unit tests for Polyderon construction & operator +=
BOOST_AUTO_TEST_CASE(PolyhedronTest) {
  std::vector<Vector3D> tvertices = {Vector3D(-1, -1, 0.), Vector3D(1., -1, 0.),
                                     Vector3D(0., 1., 0.)};
  std::vector<std::vector<size_t>> tfaces = {{0, 1, 2}};

  Polyhedron triangle(tvertices, tfaces, tfaces);
  BOOST_CHECK(tvertices == triangle.vertices);
  BOOST_CHECK(tfaces == triangle.faces);
  BOOST_CHECK(tfaces == triangle.triangularMesh);

  ObjVisualization objVis;
  GeometryView::drawPolyhedron(objVis, triangle);
  objVis.write("Polyhedron_Triangle");
  objVis.clear();

  std::vector<Vector3D> rvertices = {Vector3D(-1, -2, 0.), Vector3D(1., -2, 0.),
                                     Vector3D(1., -1., 0.),
                                     Vector3D(-1., -1., 0.)};
  std::vector<std::vector<size_t>> rfaces = {{0, 1, 2, 3}};
  std::vector<std::vector<size_t>> rmesh = {{0, 1, 2}, {2, 3, 0}};
  Polyhedron rectangle(rvertices, rfaces, rmesh);
  BOOST_CHECK(rvertices == rectangle.vertices);
  BOOST_CHECK(rfaces == rectangle.faces);
  BOOST_CHECK(rmesh == rectangle.triangularMesh);

  GeometryView::drawPolyhedron(objVis, rectangle);
  objVis.write("Polyhedron_Rectangle");
  objVis.clear();

  // Now add them
  Polyhedron tr;
  tr.merge(triangle);
  BOOST_CHECK(tr.vertices == triangle.vertices);
  BOOST_CHECK(tr.faces == triangle.faces);
  BOOST_CHECK(tr.triangularMesh == triangle.triangularMesh);
  tr.merge(rectangle);

  GeometryView::drawPolyhedron(objVis, tr);
  objVis.write("Polyhedron_TriangleRectangle");
  objVis.clear();
}

/// Unit tests for Polyderon construction & operator +=
BOOST_AUTO_TEST_CASE(PolyhedronExtent) {
  // Test a rectangle in x-y plane (at z == 0)
  std::vector<Vector3D> rvertices = {Vector3D(-1, -2, 0.), Vector3D(1., -2, 0.),
                                     Vector3D(1., -1., 0.),
                                     Vector3D(-1., -1., 0.)};

  std::vector<std::vector<size_t>> rfaces = {{0, 1, 2, 3}};
  std::vector<std::vector<size_t>> rmesh = {{0, 1, 2}, {2, 3, 0}};
  Polyhedron rectangle(rvertices, rfaces, rmesh);

  auto rExtent = rectangle.extent();
  CHECK_CLOSE_ABS(rExtent.min(binX), -1., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binX), 1., 1e-6);
  CHECK_CLOSE_ABS(rExtent.min(binY), -2., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binY), -1., 1e-6);
  CHECK_CLOSE_ABS(rExtent.min(binZ), 0., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binZ), 0., 1e-6);
  CHECK_CLOSE_ABS(rExtent.min(binR), 1., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binR), VectorHelpers::perp(rvertices[0]), 1e-6);
  CHECK_CLOSE_ABS(rExtent.min(binPhi), VectorHelpers::phi(rvertices[3]), 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binPhi), VectorHelpers::phi(rvertices[2]), 1e-6);

  // Now shift the Extent
  Vector3D shift(-1., 0., 1.);
  Transform3D shiftedTransform = Transform3D::Identity();
  shiftedTransform.pretranslate(shift);
  rExtent = rectangle.extent(shiftedTransform);
  CHECK_CLOSE_ABS(rExtent.min(binX), -2., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binX), 0., 1e-6);
  CHECK_CLOSE_ABS(rExtent.min(binY), -2., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binY), -1., 1e-6);
  CHECK_CLOSE_ABS(rExtent.min(binZ), 1., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binZ), 1., 1e-6);

  // Test a rectangle in yz - pane (at x == 3)
  rvertices = {Vector3D(3_mm, -5_mm, -10_mm), Vector3D(3_mm, 5_mm, -10_mm),
               Vector3D(3_mm, 5_mm, 10_mm), Vector3D(3_mm, -5_mm, 10_mm)};

  rectangle = Polyhedron(rvertices, rfaces, rmesh);
  rExtent = rectangle.extent();
  CHECK_CLOSE_ABS(rExtent.min(binX), 3., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binX), 3., 1e-6);
  CHECK_CLOSE_ABS(rExtent.min(binY), -5., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binY), 5., 1e-6);
  CHECK_CLOSE_ABS(rExtent.min(binZ), -10., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binZ), 10., 1e-6);
  CHECK_CLOSE_ABS(rExtent.min(binR), 3., 1e-6);
  CHECK_CLOSE_ABS(rExtent.max(binR), std::sqrt(9. + 25.), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts