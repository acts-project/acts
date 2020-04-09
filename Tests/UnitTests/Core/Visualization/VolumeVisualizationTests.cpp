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
  VolumeVisualization::test(obj, false, "");
  VolumeVisualization::test(obj, true, "_3Mesh");
}

BOOST_AUTO_TEST_CASE(VolumeVisualizationPly) {
  PlyVisualization ply;
  VolumeVisualization::test(ply, false, "");
  VolumeVisualization::test(ply, true, "_3Mesh");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts