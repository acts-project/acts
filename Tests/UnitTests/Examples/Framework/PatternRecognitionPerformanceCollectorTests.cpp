// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/Validation/PatternRecognitionPerformanceCollector.hpp"

#include <memory>
#include <numbers>

using namespace Acts;
using namespace ActsExamples;

BOOST_AUTO_TEST_CASE(InclusiveFakeRatio) {
  auto trackBackend = std::make_shared<VectorTrackContainer>();
  auto stateBackend = std::make_shared<VectorMultiTrajectory>();
  ActsExamples::TrackContainer tracks(trackBackend, stateBackend);
  TrackParticleMatching matching;
  const auto surface = Surface::makeShared<PerigeeSurface>(Vector3::Zero());

  for (int i = 0; i < 6; ++i) {
    auto track = tracks.makeTrack();
    track.parameters().setZero();
    track.parameters()[eBoundPhi] = i == 3 ? 1. : 0.;
    track.parameters()[eBoundTheta] = std::numbers::pi / 2;
    track.parameters()[eBoundQOverP] = 1.;
    track.covariance().setIdentity();
    if (i < 4) {
      track.setReferenceSurface(surface);
    }
    if (i == 0) {
      matching[track.index()].classification =
          TrackMatchClassification::Matched;
    } else if (i == 1) {
      matching[track.index()].classification =
          TrackMatchClassification::Duplicate;
    } else if (i == 2 || i == 4) {
      matching[track.index()].classification = TrackMatchClassification::Fake;
    }
  }

  ConstTrackContainer constTracks(
      std::make_shared<ConstVectorTrackContainer>(*trackBackend),
      std::make_shared<ConstVectorMultiTrajectory>(*stateBackend));
  PatternRecognitionPerformanceCollector collector(
      {}, getDefaultLogger("Test", Logging::WARNING));
  collector.fill(GeometryContext::dangerouslyDefaultConstruct(), constTracks,
                 {}, matching, {}, {});

  for (const auto& [name, efficiency] :
       collector.fakePlotTool().efficiencies()) {
    BOOST_TEST_CONTEXT(name) {
      double accepted = 0;
      double total = 0;
      for (const auto& bin : efficiency.acceptedHistogram()) {
        accepted += bin.value();
      }
      for (const auto& bin : efficiency.totalHistogram()) {
        total += bin.value();
      }
      BOOST_CHECK_EQUAL(accepted, 1.);
      BOOST_CHECK_EQUAL(total, 4.);
    }
  }
  const auto& phi =
      collector.fakePlotTool().efficiencies().at("fakeRatio_vs_phi");
  const auto unassociatedBin = phi.totalHistogram().axis(0).index(1.);
  BOOST_CHECK_EQUAL(phi.totalHistogram().at(unassociatedBin).value(), 1.);
  BOOST_CHECK_EQUAL(phi.acceptedHistogram().at(unassociatedBin).value(), 0.);
  const auto& stats = collector.stats();
  BOOST_CHECK_EQUAL(stats.nTotalTracks, 6u);
  BOOST_CHECK_EQUAL(stats.nTotalFakeTracks, 2u);
  BOOST_CHECK_EQUAL(stats.nTotalDuplicateTracks, 1u);
  BOOST_CHECK_EQUAL(stats.nTotalUnmatchedTracks, 2u);
  BOOST_CHECK_EQUAL(stats.nTotalTracksMissingRefSurface, 2u);
}
