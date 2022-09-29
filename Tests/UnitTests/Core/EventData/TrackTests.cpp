// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/Track.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

namespace {

using namespace Acts;
using namespace Acts::HashedStringLiteral;
// using namespace Acts::UnitLiterals;
// using namespace Acts::Test;
using MultiTrajectoryTraits::IndexType;
namespace bd = boost::unit_test::data;
}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataTrack)

BOOST_AUTO_TEST_CASE(Build) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto t = tc.getTrack(tc.addTrack());
  t.component<IndexType, "tipIndex"_hash>() = 5;

  BOOST_CHECK_EQUAL((t.component<IndexType, "tipIndex"_hash>()), 5);
}

BOOST_AUTO_TEST_SUITE_END()
