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

#include "Acts/Experimental/CylindricalContainerHelper.hpp"

#include <memory>
#include <vector>

#include "GeometryHelper.hpp"

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(DetectorVolumesInR) {
  auto volumesInR = createCentralDetector();
  BOOST_CHECK(volumesInR != nullptr);

  // Let's check the volume portal consistency
  const auto& portals = volumesInR->portals();
  BOOST_CHECK(portals.size() == 3u);

  const auto& volumes = volumesInR->volumes();
  BOOST_CHECK(volumes.size() == 4u);

  GeometryContext gctx = GeometryContext();

  // Get the beam pipe
  auto beamPipe = volumesInR->lowest(gctx, Vector3(0., 0., 0.));
  auto layer0 = volumesInR->lowest(gctx, Vector3(30., 0., 0.));
  auto gap = volumesInR->lowest(gctx, Vector3(50., 0., 0.));
  auto layer1 = volumesInR->lowest(gctx, Vector3(70., 0., 0.));

  // They should all be different
  BOOST_CHECK(beamPipe != layer0);
  BOOST_CHECK(layer0 != gap);
  BOOST_CHECK(gap != layer1);
  BOOST_CHECK(layer0 != layer1);

  // Let us cross-check the portals
  // -- all negative portals should be identicial and
  //    be identical with the container portal
  BOOST_CHECK(volumesInR->portals()[0] == beamPipe->portals()[0]);
  BOOST_CHECK(volumesInR->portals()[0] == layer0->portals()[0]);
  BOOST_CHECK(volumesInR->portals()[0] == gap->portals()[0]);
  BOOST_CHECK(volumesInR->portals()[0] == layer1->portals()[0]);
  // -- all positive portals should be identicial and
  //    be identical with the container portal
  BOOST_CHECK(volumesInR->portals()[1] == beamPipe->portals()[1]);
  BOOST_CHECK(volumesInR->portals()[1] == layer0->portals()[1]);
  BOOST_CHECK(volumesInR->portals()[1] == gap->portals()[1]);
  BOOST_CHECK(volumesInR->portals()[1] == layer1->portals()[1]);
  // -- the container outer portal should be the outer portal
  //    of volume 'layer1
  BOOST_CHECK(volumesInR->portals()[2] == layer1->portals()[2]);
}

BOOST_AUTO_TEST_CASE(DetectorVolumesInZ) {
  // auto volumesInR = createEndcapDetector();
  // BOOST_CHECK(volumesInR != nullptr);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
