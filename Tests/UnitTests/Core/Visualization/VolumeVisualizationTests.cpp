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
#include "VolumeVisualizationBase.hpp"

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Visualization)

BOOST_AUTO_TEST_CASE(VolumeVisualizationObj) {
  ObjVisualization obj;
  auto objCCount = VolumeVisualization::test(obj, false, "");
  BOOST_TEST(objCCount == 31877);

  auto obj3MCCount = VolumeVisualization::test(obj, true, "_3Mesh");
  BOOST_TEST(obj3MCCount == obj3MCCount);
}

BOOST_AUTO_TEST_CASE(VolumeVisualizationPly) {
  PlyVisualization ply;
  auto plyCCount = VolumeVisualization::test(ply, false, "");
  BOOST_TEST(plyCCount == 38520);

  auto ply3MCCount = VolumeVisualization::test(ply, true, "_3Mesh");
  BOOST_TEST(ply3MCCount == 38520);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts