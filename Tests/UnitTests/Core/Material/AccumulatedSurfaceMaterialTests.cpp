// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/AccumulatedSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"

#include <climits>

namespace Acts {
namespace Test {

/// Test the constructors
BOOST_AUTO_TEST_CASE(AccumulatedSurfaceMaterial_construction_test) {
  // Test:
  // HomogeneousSurfaceMaterial accumulation
  AccumulatedSurfaceMaterial material0D{};
  auto accMat0D = material0D.accumulatedMaterial();
  BOOST_CHECK_EQUAL(accMat0D.size(), 1u);
  BOOST_CHECK_EQUAL(accMat0D[0].size(), 1u);
  BOOST_CHECK_EQUAL(material0D.splitFactor(), 0.);

  // Test:
  // BinnesSurfaceMatieral accumulation - 1D
  BinUtility binUtility1D(10, -5., 5., open, binX);
  AccumulatedSurfaceMaterial material1D{binUtility1D};
  auto accMat1D = material1D.accumulatedMaterial();
  BOOST_CHECK_EQUAL(accMat1D.size(), 1u);
  BOOST_CHECK_EQUAL(accMat1D[0].size(), 10u);

  // Test:
  // BinnesSurfaceMatieral accumulation - 2D
  BinUtility binUtility2D(10, -5., 5., open, binX);
  binUtility2D += BinUtility(20, -10., 10., open, binY);
  AccumulatedSurfaceMaterial material2D{binUtility2D};
  auto accMat2D = material2D.accumulatedMaterial();
  BOOST_CHECK_EQUAL(accMat2D.size(), 20u);
  for (size_t ib = 0; ib < accMat2D.size(); ++ib) {
    BOOST_CHECK_EQUAL(accMat2D[ib].size(), 10u);
  }
}

/// Test the filling and conversion
BOOST_AUTO_TEST_CASE(AccumulatedSurfaceMaterial_fill_convert_0D) {
  Material mat = Material::fromMolarDensity(1., 1., 1., 1., 1.);
  MaterialSlab one(mat, 1.);
  MaterialSlab two(mat, 2.);

  AccumulatedSurfaceMaterial material0D{};
  // assign 2 one steps
  material0D.accumulate(Vector2D{0., 0.}, one);
  material0D.accumulate(Vector2D{0., 0.}, one);
  material0D.trackAverage();
  // assign 1 double step
  material0D.accumulate(Vector3D(0., 0., 0.), two);
  material0D.trackAverage();
  // get the single matrix
  auto accMat0D = material0D.accumulatedMaterial();
  auto accMatProp0D = accMat0D[0][0];
  auto [matProp0D, trackCount] = accMatProp0D.totalAverage();

  BOOST_CHECK_EQUAL(matProp0D.thicknessInX0(), two.thicknessInX0());
  BOOST_CHECK_EQUAL(trackCount, 2u);
}

/// Test the filling and conversion
BOOST_AUTO_TEST_CASE(AccumulatedSurfaceMaterial_fill_convert_1D) {
  Material mat = Material::fromMolarDensity(1., 1., 1., 1., 1.);
  MaterialSlab one(mat, 1.);
  MaterialSlab two(mat, 2.);
  MaterialSlab three(mat, 3.);
  MaterialSlab four(mat, 4.);

  // BinnesSurfaceMatieral accumulation - 2D
  BinUtility binUtility2D(2, -1., 1., open, binX);
  binUtility2D += BinUtility(2, -1., 1., open, binY);
  AccumulatedSurfaceMaterial material2D{binUtility2D};

  // assign in the different bins
  // event 0
  material2D.accumulate(Vector2D{-0.5, -0.5}, one);
  material2D.accumulate(Vector2D{0.5, -0.5}, two);
  material2D.accumulate(Vector2D{-0.5, 0.5}, three);
  material2D.accumulate(Vector2D{0.5, 0.5}, four);
  material2D.trackAverage();
  // event 1
  material2D.accumulate(Vector2D{0.5, -0.5}, two);
  material2D.accumulate(Vector2D{-0.5, 0.5}, three);
  material2D.accumulate(Vector2D{0.5, 0.5}, four);
  material2D.trackAverage();
  // event 2
  material2D.accumulate(Vector2D{-0.5, 0.5}, three);
  material2D.accumulate(Vector2D{0.5, 0.5}, four);
  material2D.trackAverage();
  // event 2
  material2D.accumulate(Vector2D{0.5, 0.5}, four);
  material2D.trackAverage();
  // get the single matrix
  auto accMat2D = material2D.accumulatedMaterial();
  // the accumulated properties
  auto [accMatProp00, trackCount00] = accMat2D[0][0].totalAverage();
  auto [accMatProp01, trackCount01] = accMat2D[0][1].totalAverage();
  auto [accMatProp10, trackCount10] = accMat2D[1][0].totalAverage();
  auto [accMatProp11, trackCount11] = accMat2D[1][1].totalAverage();

  BOOST_CHECK_EQUAL(accMatProp00.thicknessInX0(), one.thicknessInX0());
  BOOST_CHECK_EQUAL(accMatProp01.thicknessInX0(), two.thicknessInX0());
  BOOST_CHECK_EQUAL(accMatProp10.thicknessInX0(), three.thicknessInX0());
  BOOST_CHECK_EQUAL(accMatProp11.thicknessInX0(), four.thicknessInX0());
  BOOST_CHECK_EQUAL(trackCount00, 1u);
  BOOST_CHECK_EQUAL(trackCount01, 2u);
  BOOST_CHECK_EQUAL(trackCount10, 3u);
  BOOST_CHECK_EQUAL(trackCount11, 4u);
}

}  // namespace Test
}  // namespace Acts
