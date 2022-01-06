// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Experimental/InternalBlueprint.hpp"
#include "Acts/Experimental/SurfaceLinks.hpp"

#include "GeometryHelper.hpp"

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Experimental)

// Single Link test
BOOST_AUTO_TEST_CASE(InternalBlueprintConstruction) {
  // These are two cylinders
  auto layer0 = surfacesCylinder(8.4, 36., 0.145, 32., 2., 5., {16, 14});
  auto layer1 = surfacesCylinder(8.4, 36., 0.145, 72., 2., 5., {32, 14});

  GeometryContext gctx;

  GeometricExtent layer0Extent;
  GeometricExtent layer1Extent;

  // Create the two Blue prints for the layers
  InternalBlueprint layer0Thight(gctx, layer0, AllInternalSurfaces{},
                                 {{binZ, {0., 0.}}, {binR, {0., 0.}}},
                                 "layer0");
  InternalBlueprint layer1Thight(gctx, layer1, AllInternalSurfaces{},
                                 {{binZ, {0., 0.}}, {binR, {0., 0.}}},
                                 "layer1");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts