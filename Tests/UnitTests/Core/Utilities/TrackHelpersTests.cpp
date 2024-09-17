// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Utilities)

BOOST_AUTO_TEST_CASE(CalculateQuantities) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto t = tc.makeTrack();

  auto ts = t.appendTrackState();
  ts.typeFlags().set(MeasurementFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(OutlierFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(MeasurementFlag);
  ts.typeFlags().set(SharedHitFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(HoleFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(OutlierFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(HoleFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(MeasurementFlag);
  ts.typeFlags().set(SharedHitFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(OutlierFlag);

  calculateTrackQuantities(t);
  BOOST_CHECK_EQUAL(t.nHoles(), 2);
  BOOST_CHECK_EQUAL(t.nMeasurements(), 3);
  BOOST_CHECK_EQUAL(t.nOutliers(), 3);
  BOOST_CHECK_EQUAL(t.nSharedHits(), 2);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
