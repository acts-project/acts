// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/PortalNavigation.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <memory>
#include <stdexcept>
#include <vector>

namespace Acts::Experimental {
class DetectorVolume {};
}  // namespace Acts::Experimental

using namespace Acts;

// A test context
GeometryContext tContext;

auto volumeA = std::make_shared<Experimental::DetectorVolume>();
auto volumeB = std::make_shared<Experimental::DetectorVolume>();
auto volumeC = std::make_shared<Experimental::DetectorVolume>();
auto volumeD = std::make_shared<Experimental::DetectorVolume>();

Experimental::NavigationState nState;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(NavigationSuite)

// These tests check the behavior of the volume updators, i.e. the
// helper delegates that set/reset the volume raw pointer in the
// NavigaitonState according to some given information.
//
BOOST_AUTO_TEST_CASE(UnconnectedUpdate) {
  Experimental::ExternalNavigationDelegate ucUpdater;
  BOOST_CHECK(!ucUpdater.connected());
}

// The end of world is reached
BOOST_AUTO_TEST_CASE(EndOfWorldUpdate) {
  nState.currentVolume = volumeA.get();
  BOOST_CHECK_EQUAL(nState.currentVolume, volumeA.get());

  Experimental::EndOfWorld eow;
  eow.update(tContext, nState);

  BOOST_CHECK_EQUAL(nState.currentVolume, nullptr);
}

// A single link exists and this is set
BOOST_AUTO_TEST_CASE(SingleVolumeUpdate) {
  nState.currentVolume = volumeA.get();
  BOOST_CHECK_EQUAL(nState.currentVolume, volumeA.get());

  Experimental::SingleDetectorVolumeNavigation svu(volumeB.get());
  svu.update(tContext, nState);

  BOOST_CHECK_EQUAL(nState.currentVolume, volumeB.get());

  BOOST_CHECK_THROW(Experimental::SingleDetectorVolumeNavigation(nullptr),
                    std::invalid_argument);
}

// A typlical volume array in 1 dimension (bound, not closed)
BOOST_AUTO_TEST_CASE(VolumeArrayUpdate) {
  std::vector<double> zArray = {-200, -100, 100, 400, 1000};

  std::vector<const Experimental::DetectorVolume*> volumes = {
      volumeA.get(), volumeB.get(), volumeC.get(), volumeD.get()};
  Experimental::BoundVolumesGrid1Navigation bvg(zArray, AxisDirection::AxisZ,
                                                volumes);
  // Reset the navigation state
  nState.currentVolume = nullptr;

  // Check the volume retrieval
  nState.position = Vector3(0., 0., -150.);
  bvg.update(tContext, nState);
  BOOST_CHECK_EQUAL(nState.currentVolume, volumeA.get());

  nState.position = Vector3(0., 0., 600.);
  bvg.update(tContext, nState);
  BOOST_CHECK_EQUAL(nState.currentVolume, volumeD.get());

  // Check a shifted one
  Transform3 shift300 = Transform3::Identity();
  shift300.pretranslate(Vector3(0, 0, 300));

  Experimental::BoundVolumesGrid1Navigation bvgs(zArray, AxisDirection::AxisZ,
                                                 volumes, shift300.inverse());

  // 150 (-300) -> transforms to -150, hence it yields A
  nState.position = Vector3(0., 0., 150.);
  bvgs.update(tContext, nState);
  BOOST_CHECK_EQUAL(nState.currentVolume, volumeA.get());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
