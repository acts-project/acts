// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <vector>

namespace Acts::Experimental {
/// Define a dummy detector volume
class DetectorVolume {};
}  // namespace Acts::Experimental

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(NavigationSuite)

BOOST_AUTO_TEST_CASE(NavigationState) {
  // A navigation state struct
  Experimental::NavigationState nState;
  auto dTransform = Transform3::Identity();

  // A rectangle bound surface
  auto rectangle = std::make_shared<RectangleBounds>(10., 100.);
  auto surfaceA = Surface::makeShared<PlaneSurface>(dTransform, rectangle);
  auto surfaceB = Surface::makeShared<PlaneSurface>(dTransform, rectangle);
  auto surfaceC = Surface::makeShared<PlaneSurface>(dTransform, rectangle);

  // portal surfaces
  auto pSurfaceA = Surface::makeShared<PlaneSurface>(dTransform, rectangle);
  auto pSurfaceB = Surface::makeShared<PlaneSurface>(dTransform, rectangle);

  // Create a few fake portals out of it
  auto portalA = std::make_shared<Experimental::Portal>(pSurfaceA);
  auto portalB = std::make_shared<Experimental::Portal>(pSurfaceB);

  std::vector<const Surface*> surfaces = {surfaceA.get(), surfaceB.get(),
                                          surfaceC.get()};
  std::vector<const Experimental::Portal*> portals = {portalA.get(),
                                                      portalB.get()};

  auto dVolume = std::make_unique<Experimental::DetectorVolume>();
  const auto volume = dVolume.get();

  Experimental::DetectorVolumeFiller::fill(nState, volume);
  BOOST_CHECK_EQUAL(nState.currentVolume, volume);

  Experimental::PortalsFiller::fill(nState, portals);
  BOOST_CHECK_EQUAL(nState.surfaceCandidates.size(), portals.size());

  Experimental::SurfacesFiller::fill(nState, surfaces);
  BOOST_CHECK_EQUAL(nState.surfaceCandidates.size(),
                    portals.size() + surfaces.size());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
