// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>

#include "Acts/Visualization/IVisualization.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"
#include "Acts/Visualization/PlyVisualization.hpp"
#include "VisualizationTester.hpp"
#include "VolumeViewBase.hpp"

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Visualization)

// Tests if the written obj output is well formatted
BOOST_AUTO_TEST_CASE(VolumeViewObj) {
  ObjVisualization obj;

  // Standard test
  bool triangulate = false;
  auto objTest = VolumeViewTest::run(obj, triangulate, "");
  auto objErrors = testObjString(objTest, triangulate);
  std::cout << "Volumes Obj Test    : " << objTest.size()
            << " characters written with " << objErrors.size() << " errors."
            << std::endl;
  BOOST_CHECK(objErrors.size() == 0);
  for (auto objerr : objErrors) {
    std::cout << objerr << std::endl;
  }

  // Triangular mesh test
  triangulate = true;
  auto objTest3M = VolumeViewTest::run(obj, triangulate, "_3M");
  auto objErrors3M = testObjString(objTest3M, triangulate);
  std::cout << "Volumes Obj Test 3M : " << objTest3M.size()
            << " characters written with " << objErrors3M.size() << " errors."
            << std::endl;
  BOOST_CHECK(objErrors3M.size() == 0);
  for (auto objerr : objErrors3M) {
    std::cout << objerr << std::endl;
  }
}

// Tests if the written ply output is well formatted
BOOST_AUTO_TEST_CASE(VolumeViewPly) {
  PlyVisualization ply;

  // Standard test
  bool triangulate = false;
  auto plyTest = VolumeViewTest::run(ply, triangulate, "");
  auto plyErrors = testPlyString(plyTest, triangulate);
  std::cout << "Volumes Ply Test    : " << plyTest.size()
            << " characters written with " << plyErrors.size() << " errors."
            << std::endl;
  BOOST_CHECK(plyErrors.size() == 0);
  for (auto plyerr : plyErrors) {
    std::cout << plyerr << std::endl;
  }

  // Triangular mesh test
  triangulate = true;
  auto plyTest3M = VolumeViewTest::run(ply, triangulate, "_3M");
  auto plyErrors3M = testPlyString(plyTest3M, triangulate);
  std::cout << "Volumes Ply Test 3M : " << plyTest3M.size()
            << " characters written with " << plyErrors3M.size() << " errors."
            << std::endl;
  BOOST_CHECK(plyErrors3M.size() == 0);
  for (auto plyerr : plyErrors3M) {
    std::cout << plyerr << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts