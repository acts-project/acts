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

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/PortalLinks.hpp"
#include "Acts/Experimental/VolumeLinks.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Utilities/detail/Axis.hpp"

#include <memory>
#include <vector>

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Experimental)

// Single Link test
BOOST_AUTO_TEST_CASE(SingleAndMultiplePortalLink) {
  // Create some unique bounds
  auto cylVolBounds0 =
      std::make_unique<Acts::CylinderVolumeBounds>(1., 4., 10.);
  auto volume0 = DetectorVolume::makeShared(
      Transform3::Identity(), std::move(cylVolBounds0), "CylinderVolume0");

  SinglePortalLink spLink0{volume0.get()};
  BOOST_CHECK(spLink0.volume == volume0.get());

  auto cylVolBounds1 =
      std::make_unique<Acts::CylinderVolumeBounds>(4., 6., 10.);
  auto volume1 = DetectorVolume::makeShared(
      Transform3::Identity(), std::move(cylVolBounds1), "CylinderVolume1");
  SinglePortalLink spLink1{volume1.get()};
  BOOST_CHECK(spLink1.volume == volume1.get());

  std::vector<ActsScalar> rBoundariesInner = {1., 4., 6.};
  VariableVolumeLink rInnerLink(detail::VariableAxis(rBoundariesInner), binR);

  // Let us test the mulit portal link
  Vector3 getFirst(2., 0., 0.);
  Vector3 getSecond(5., 0., 0.);

  unsigned int i0 = rInnerLink(Transform3::Identity(), getFirst);
  BOOST_CHECK(i0 == 0u);

  unsigned int i1 = rInnerLink(Transform3::Identity(), getSecond);
  BOOST_CHECK(i1 == 1u);

  MultiplePortalLink rInnerMultiPortalLink{{spLink0, spLink1},
                                           std::move(rInnerLink)};
  Vector3 atFirst(1., 0., 0.);
  Vector3 atSecond(1., 0., 0.);
  Vector3 direction(1., 0., 0.);
  GeometryContext gctx = GeometryContext();

  const auto& portals = volume1->portals();
  const Portal& firstPortal = (*portals[3]);

  auto environment =
      rInnerMultiPortalLink(gctx, firstPortal, atFirst, direction, true);
  BOOST_CHECK(environment.volume == volume0.get());
}

// This creates Single links by construction
BOOST_AUTO_TEST_CASE(DetectorVolumePortalCandidates) {
  GeometryContext gctx;

  // Create some unique bounds
  auto cylVolBounds = std::make_unique<Acts::CylinderVolumeBounds>(1., 4., 10.);
  auto volume = DetectorVolume::makeShared(
      Transform3::Identity(), std::move(cylVolBounds), "CylinderVolume");

  // Get the portal raw pointers
  auto portals = volume->portals();

  // Tst the portals from an entry point
  Vector3 atInner(1., 0., 0.);
  Vector3 direction = Vector3(1., 1., 1.).normalized();

  auto pCandidates = portalCandidates(gctx, volume->portals(), atInner, direction);

  BOOST_TEST(pCandidates.size(), 4u);
  // Exit should be through the outer tube
  BOOST_TEST(pCandidates[0].object, portals[2]);
  // Next candidate is positize z
  BOOST_TEST(pCandidates[1].object, portals[1]);
  // Then negative z
  BOOST_TEST(pCandidates[2].object, portals[0]);
  // Finally, on surface is ranked last
  BOOST_TEST(pCandidates[3].object, portals[3]);

  // Test the portals from an inside point (same setup)
  Vector3 atEpsilonFromInner(1.5, 0., 0.);
  // Result should be unchanged
  pCandidates =
      portalCandidates(gctx, volume->portals(), atEpsilonFromInner, direction);

  BOOST_TEST(pCandidates.size(), 4u);
  // Exit should be through the outer tube
  BOOST_TEST(pCandidates[0].object, portals[2]);
  // Next candidate is positize z
  BOOST_TEST(pCandidates[1].object, portals[1]);
  // Then negative z
  BOOST_TEST(pCandidates[2].object, portals[0]);
  // Finally, on surface is ranked last
  BOOST_TEST(pCandidates[3].object, portals[3]);

  // Test from inner, but reverse
  pCandidates = portalCandidates(gctx, volume->portals(), atEpsilonFromInner,
                                 -1. * direction);

  BOOST_TEST(pCandidates.size(), 4u);
  // Exit should be through the inner tube
  BOOST_TEST(pCandidates[0].object, portals[3]);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
