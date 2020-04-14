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
#include "SurfaceVisualizationBase.hpp"

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Visualization)

/// The tests in this section are regression tests only in order
/// to catch any unexpected changes in the output format.
///
BOOST_AUTO_TEST_CASE(SurfaceVisualizationObj) {
  ObjVisualization obj;
  size_t objCCount = SurfaceVisualization::test(obj, false, "");
  // Due to rounding errors, there can be a -0 or 0 printed
  // allow small tolerance for the moment
  BOOST_CHECK_CLOSE(1. * objCCount, 1. * 56360, 0.1);

  size_t objC3MCount = SurfaceVisualization::test(obj, true, "_3Mesh");
  BOOST_CHECK_EQUAL(objC3MCount, 66878);
}

/// The tests in this section are regression tests only in order
/// to catch any unexpected changes in the output format.
///
BOOST_AUTO_TEST_CASE(SurfaceVisualizationPly) {
  PlyVisualization ply;
  size_t plyCCount = SurfaceVisualization::test(ply, false, "");
  BOOST_CHECK_EQUAL(plyCCount, 76912);

  size_t plyC3MCount = SurfaceVisualization::test(ply, true, "_3Mesh");
  BOOST_CHECK_EQUAL(plyC3MCount, 76912);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts