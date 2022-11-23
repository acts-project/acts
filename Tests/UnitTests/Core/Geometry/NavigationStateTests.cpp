// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationState.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>

// A test context
Acts::GeometryContext tContext;

namespace Acts {
namespace Experimental {
/// Define a dummy detector volume
class DetectorVolume {};
}  // namespace Experimental
}  // namespace Acts

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(NavigationState) {
  // A navigation state struct
  Acts::Experimental::NavigationState nState;
  auto dTransform = Acts::Transform3::Identity();

  // A rectangle bound surface
  auto rectangle = std::make_shared<Acts::RectangleBounds>(10., 100.);
  auto surfaceA =
      Acts::Surface::makeShared<Acts::PlaneSurface>(dTransform, rectangle);
  auto surfaceB =
      Acts::Surface::makeShared<Acts::PlaneSurface>(dTransform, rectangle);
  auto surfaceC =
      Acts::Surface::makeShared<Acts::PlaneSurface>(dTransform, rectangle);

  // portal surfaces
  auto pSurfaceA =
      Acts::Surface::makeShared<Acts::PlaneSurface>(dTransform, rectangle);
  auto pSurfaceB =
      Acts::Surface::makeShared<Acts::PlaneSurface>(dTransform, rectangle);

  // Create a few fake portals out of it
  auto portalA = Acts::Experimental::Portal::makeShared(pSurfaceA);
  auto portalB = Acts::Experimental::Portal::makeShared(pSurfaceB);

  std::vector<const Acts::Surface*> surfaces = {surfaceA.get(), surfaceB.get(),
                                                surfaceC.get()};
  std::vector<const Acts::Experimental::Portal*> portals = {portalA.get(),
                                                            portalB.get()};

  auto dVolume = std::make_unique<Acts::Experimental::DetectorVolume>();
  const auto volume = dVolume.get();

  Acts::Experimental::DetectorVolumeFiller::fill(nState, volume);
  BOOST_CHECK(nState.currentVolume == volume);

  Acts::Experimental::PortalsFiller::fill(nState, portals);
  BOOST_CHECK(nState.surfaceCandidates.size() == portals.size());

  Acts::Experimental::SurfacesFiller::fill(nState, surfaces);
  BOOST_CHECK(nState.surfaceCandidates.size() ==
              portals.size() + surfaces.size());
}

BOOST_AUTO_TEST_SUITE_END()
