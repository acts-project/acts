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
#include "Acts/Geometry/Extent.hpp"

#include "GeometryHelper.hpp"

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Experimental)

// Single Link test
BOOST_AUTO_TEST_CASE(InternalBlueprintConstruction) {
  auto layer0 = surfacesCylinder(8.4, 36., 0.145, 32., 2., 5., {16, 14});
  auto layer1 = surfacesCylinder(8.4, 36., 0.145, 72., 2., 5., {32, 14});

  std::vector<std::shared_ptr<Acts::Surface>> layer01;
  layer01.insert(layer01.begin(), layer0.begin(), layer0.end());
  layer01.insert(layer01.begin(), layer1.begin(), layer1.end());

  // Check that test environement is ok
  BOOST_CHECK(layer01.size() == layer0.size() + layer1.size());

  auto zeroEnvelope = std::vector<std::pair<ActsScalar, ActsScalar>>(
      size_t(binValues), {0., 0.});

  Extent layer0Ext(true);
  layer0Ext.ranges[binR] = {0., 40.};

  Extent layer1Ext(true);
  layer1Ext.ranges[binR] = {40., 100.};

  GeometryContext gctx;

  // Create the two Blue prints for the layers
  InternalBlueprint layer0Bp(gctx, layer01, AllSurfaces{}, layer0Ext, zeroEnvelope);
  InternalBlueprint layer1Bp(gctx, layer01, AllSurfaces{}, layer1Ext, zeroEnvelope);

  BOOST_CHECK(layer0Bp.surfaces().size() == layer0.size());
  BOOST_CHECK(layer1Bp.surfaces().size() == layer1.size());

}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts