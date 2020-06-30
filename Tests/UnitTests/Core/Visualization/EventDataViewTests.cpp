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
#include "EventDataViewBase.hpp"
#include "VisualizationTester.hpp"

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Visualization)

BOOST_AUTO_TEST_CASE(BoundParametersVisualizationObj) {
  ObjVisualization obj;
  auto objTest = EventDataViewTest::testBoundParameters(obj);
  auto objErrors = testObjString(objTest);
  BOOST_CHECK(objErrors.size() == 0);
  for (auto objerr : objErrors) {
    std::cout << objerr << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(BoundParametersVisualizationPly) {
  PlyVisualization ply;
  auto plyTest = EventDataViewTest::testBoundParameters(ply);
  auto plyErrors = testPlyString(plyTest);
  BOOST_CHECK(plyErrors.size() == 0);
  for (auto plyerr : plyErrors) {
    std::cout << plyerr << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryVisualizationObj) {
  ObjVisualization obj;
  auto objTest = EventDataViewTest::testMultiTrajectory(obj);
  auto objErrors = testObjString(objTest);
  BOOST_CHECK(objErrors.size() == 0);
  for (auto objerr : objErrors) {
    std::cout << objerr << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryVisualizationPly) {
  PlyVisualization ply;
  auto plyTest = EventDataViewTest::testMultiTrajectory(ply);
  auto plyErrors = testPlyString(plyTest);
  BOOST_CHECK(plyErrors.size() == 0);
  for (auto plyerr : plyErrors) {
    std::cout << plyerr << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
