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
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Enumerate.hpp"
#include "Acts/Experimental/Portal.hpp"

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

  // The outer portal of the beamPipe should be the inner of the layer0
  BOOST_CHECK(beamPipe->portals()[2] == layer0->portals()[3]);
  // And so on ..
  BOOST_CHECK(layer0->portals()[2] == gap->portals()[3]);
  BOOST_CHECK(gap->portals()[2] == layer1->portals()[3]);


  // -- the container outer portal should be the outer portal
  //    of volume 'layer1
  BOOST_CHECK(volumesInR->portals()[2] == layer1->portals()[2]);
}

BOOST_AUTO_TEST_CASE(DetectorVolumesInZ) {
  // Get the volume from the helper
  auto volumesInZ = createEndcapDetector();
  BOOST_CHECK(volumesInZ != nullptr);

  // Count the portals
  const auto& portals = volumesInZ->portals();
  BOOST_CHECK(portals.size() == 3u);

  const auto& volumes = volumesInZ->volumes();
  BOOST_CHECK(volumes.size() == 3u);

  GeometryContext gctx = GeometryContext();

  // Get the beam pipe
  auto layer0 = volumesInZ->lowest(gctx, Vector3(0., 0., 502.5));
  auto gap = volumesInZ->lowest(gctx, Vector3(0., 0., 525.));
  auto layer1 = volumesInZ->lowest(gctx, Vector3(50., 0., 557.5));

  // They should all be different
  BOOST_CHECK(layer0 != gap);
  BOOST_CHECK(gap != layer1);
  BOOST_CHECK(layer0 != layer1);

  // Let us cross-check the portals
  // -- all outer portals should be identical
  BOOST_CHECK(layer0->portals()[2] == gap->portals()[2]);
  BOOST_CHECK(gap->portals()[2] == layer1->portals()[2]);
  BOOST_CHECK(volumesInZ->portals()[2] == layer1->portals()[2]);

  // The intermediate portals should be shared
  BOOST_CHECK(layer0->portals()[1] == gap->portals()[0]);
  BOOST_CHECK(gap->portals()[1] == layer1->portals()[0]);

}

BOOST_AUTO_TEST_CASE(Detector) {


  auto detector = createDetector();
  BOOST_CHECK(detector != nullptr);
  BOOST_CHECK(detector->name() == "Detector");

  // Get the portals and the volumes
  const auto& portals = detector->portals();
  BOOST_CHECK(portals.size() == 3u);

  const auto& volumes = detector->volumes();
  BOOST_CHECK(volumes.size() == 3u);

  GeometryContext gctx = GeometryContext();

  // Negative sector
  auto layer1Neg = detector->lowest(gctx, Vector3(0., 0., -502.5));
  auto gapNeg = detector->lowest(gctx, Vector3(0., 0., -525.));
  auto layer0Neg = detector->lowest(gctx, Vector3(50., 0., -557.5));

  BOOST_CHECK(layer0Neg->portals()[2] == gapNeg->portals()[2]);
  BOOST_CHECK(gapNeg->portals()[2] == layer1Neg->portals()[2]);

  // The intermediate portals should be shared
  BOOST_CHECK(layer0Neg->portals()[1] == gapNeg->portals()[0]);
  BOOST_CHECK(gapNeg->portals()[1] == layer1Neg->portals()[0]);

  // Central sector
  auto beamPipe = detector->lowest(gctx, Vector3(0., 0., 0.));
  auto layer0 = detector->lowest(gctx, Vector3(30., 0., 0.));
  auto gap = detector->lowest(gctx, Vector3(50., 0., 0.));
  auto layer1 = detector->lowest(gctx, Vector3(70., 0., 0.));

  // Cross-check that the portals are still equal
  BOOST_CHECK(beamPipe->portals()[0] == layer0->portals()[0]);
  BOOST_CHECK(layer0->portals()[0] == gap->portals()[0]);
  BOOST_CHECK(gap->portals()[0] == layer1->portals()[0]);

  BOOST_CHECK(beamPipe->portals()[1] == layer0->portals()[1]);
  BOOST_CHECK(layer0->portals()[1] == gap->portals()[1]);
  BOOST_CHECK(gap->portals()[1] == layer1->portals()[1]);

  // The outer portal of the beamPipe should be the inner of the layer0
  BOOST_CHECK(beamPipe->portals()[2] == layer0->portals()[3]);
  // And so on ..
  BOOST_CHECK(layer0->portals()[2] == gap->portals()[3]);
  BOOST_CHECK(gap->portals()[2] == layer1->portals()[3]);

  // Positive sector
  auto layer0Pos = detector->lowest(gctx, Vector3(0., 0., 502.5));
  auto gapPos = detector->lowest(gctx, Vector3(0., 0., 525.));
  auto layer1Pos = detector->lowest(gctx, Vector3(50., 0., 557.5));

  BOOST_CHECK(layer0Pos->portals()[2] == gapPos->portals()[2]);
  BOOST_CHECK(gapPos->portals()[2] == layer1Pos->portals()[2]);

  // The intermediate portals should be shared
  BOOST_CHECK(layer0Pos->portals()[1] == gapPos->portals()[0]);
  BOOST_CHECK(gapPos->portals()[1] == layer1Pos->portals()[0]);

  // Check that the volumes are glued together
  BOOST_CHECK(layer1Neg->portals()[1] == beamPipe->portals()[0]);
  BOOST_CHECK(layer1Neg->portals()[1] == layer0->portals()[0]);
  BOOST_CHECK(layer1Neg->portals()[1] == gap->portals()[0]);
  BOOST_CHECK(layer1Neg->portals()[1] == layer1->portals()[0]);

  BOOST_CHECK(layer0Pos->portals()[0] == beamPipe->portals()[1]);
  BOOST_CHECK(layer0Pos->portals()[0] == layer0->portals()[1]);
  BOOST_CHECK(layer0Pos->portals()[0] == gap->portals()[1]);
  BOOST_CHECK(layer0Pos->portals()[0] == layer1->portals()[1]);

  std::vector<const DetectorVolume*> navigationVolumes = {
      layer1Neg, gapNeg, layer0Neg, beamPipe, layer0,
      gap,       layer1, layer0Pos, gapPos,   layer1Pos};

  for (const auto& nv : navigationVolumes) {
    BOOST_CHECK(nv != nullptr);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts