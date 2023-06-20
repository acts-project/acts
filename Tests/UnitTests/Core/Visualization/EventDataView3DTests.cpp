// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/PlyVisualization3D.hpp"

#include <iostream>
#include <string>
#include <vector>

#include "EventDataView3DBase.hpp"
#include "Visualization3DTester.hpp"

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Visualization)

BOOST_AUTO_TEST_CASE(BoundTrackParametersVisualizationObj) {
  ObjVisualization3D obj;
  auto objTest = EventDataView3DTest::testBoundTrackParameters(obj);
  auto objErrors = testObjString(objTest);
  BOOST_CHECK(objErrors.empty());
  for (const auto& objerr : objErrors) {
    std::cout << objerr << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(BoundTrackParametersVisualizationPly) {
  PlyVisualization3D ply;
  auto plyTest = EventDataView3DTest::testBoundTrackParameters(ply);
  auto plyErrors = testPlyString(plyTest);
  BOOST_CHECK(plyErrors.empty());
  for (const auto& plyerr : plyErrors) {
    std::cout << plyerr << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryVisualizationObj) {
  ObjVisualization3D obj;
  auto objTest = EventDataView3DTest::testMultiTrajectory(obj);
  auto objErrors = testObjString(objTest);
  BOOST_CHECK(objErrors.empty());
  for (const auto& objerr : objErrors) {
    std::cout << objerr << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryVisualizationPly) {
  PlyVisualization3D ply;
  auto plyTest = EventDataView3DTest::testMultiTrajectory(ply);
  auto plyErrors = testPlyString(plyTest);
  BOOST_CHECK(plyErrors.empty());
  for (const auto& plyerr : plyErrors) {
    std::cout << plyerr << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
