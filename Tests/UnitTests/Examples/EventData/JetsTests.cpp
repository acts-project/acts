// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "ActsExamples/EventData/Jets.hpp"

using namespace Acts;
using namespace ActsExamples;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(FastJetSuite)

BOOST_AUTO_TEST_CASE(TruthParticleOneJet) {
  ActsFatras::Barcode barcode;
  ActsExamples::TruthJet jet(Acts::Vector4(100, 0, 0, 100),
                             ActsExamples::JetLabel::Unknown);
  jet.setConstituents(std::vector<ActsFatras::Barcode>{barcode});
  BOOST_CHECK_EQUAL(jet.constituents().size(), 1);
  BOOST_CHECK_EQUAL(jet.constituents()[0], barcode);
}

BOOST_AUTO_TEST_CASE(SingleTrackJet) {
  Acts::TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto track = tracks.makeTrack();

  track.parameters()[Acts::eBoundLoc0] = 10.0;
  track.parameters()[Acts::eBoundLoc1] = 0.0;
  track.parameters()[Acts::eBoundTime] = 0.0;

  auto constTrack = tracks.getTrack(0);
  AnyConstTrackProxy anyConstTrack(constTrack);

  std::vector<Acts::AnyConstTrackProxy> constituents{anyConstTrack};

  ActsExamples::TrackJet jet(Acts::Vector4(100, 0, 0, 100),
                             ActsExamples::JetLabel::Unknown);
  std::vector<Acts::AnyConstTrackProxy> jetConstituents;
  jetConstituents.push_back(anyConstTrack);
  jet.setConstituents(jetConstituents);

  BOOST_CHECK_EQUAL(jet.constituents().size(), 1);
  BOOST_CHECK(jet.constituents()[0].index() == tracks.getTrack(0).index());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
