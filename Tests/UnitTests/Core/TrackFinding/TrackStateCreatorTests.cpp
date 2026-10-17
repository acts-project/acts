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
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"
#include "Acts/Utilities/Holders.hpp"

#include <algorithm>
#include <array>
#include <numeric>

namespace {
using Container = Acts::TrackContainer<Acts::VectorTrackContainer,
                                       Acts::VectorMultiTrajectory,
                                       Acts::detail::ValueHolder>;
using Creator =
    Acts::TrackStateCreator<std::vector<Acts::SourceLink>::const_iterator,
                            Container>;
struct Calibrator {
  bool extra = false;
  void calibrate([[maybe_unused]] const Acts::GeometryContext& gctx,
                 [[maybe_unused]] const Acts::CalibrationContext& cctx,
                 const Acts::SourceLink& sl,
                 Container::TrackStateProxy state) const {
    int id = sl.get<int>();
    state.setUncalibratedSourceLink(Acts::SourceLink{sl});
    const Acts::Vector2 value{id * 0.1, id * 0.2};
    const Acts::SquareMatrix2 covariance =
        Acts::SquareMatrix2::Identity() * (id + 1);
    state.allocateCalibrated(value, covariance);
    state.setProjectorSubspaceIndices(std::array<std::uint8_t, 2>{0, 1});
    state.component<int>("marker") = id + 100;
    if (extra) {
      state.addComponents(Acts::TrackStatePropMask::Filtered);
      state.filtered().setZero();
      state.filteredCovariance().setIdentity();
    }
  }
};
struct Selector {
  std::size_t selected = 2;
  bool outlier = false;
  auto select(Creator::candidate_container_t& candidates, bool& isOutlier,
              [[maybe_unused]] const Acts::Logger& logger) const {
    std::ranges::reverse(candidates);
    isOutlier = outlier;
    using Range = std::pair<Creator::candidate_container_t::iterator,
                            Creator::candidate_container_t::iterator>;
    return Acts::Result<Range>::success(
        {candidates.begin(), candidates.begin() + selected});
  }
};
}  // namespace

BOOST_AUTO_TEST_CASE(SelectedStatesMatchOriginalCopyPath) {
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  const Acts::CalibrationContext cctx;
  auto surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  Acts::BoundVector params;
  params << 0., 0., 0.2, 1., 0.001, 0.;
  Creator::BoundState bound{
      Acts::BoundTrackParameters(surface, params, Acts::BoundMatrix::Identity(),
                                 Acts::ParticleHypothesis::pion()),
      Acts::BoundMatrix::Identity() * 2., 42.};
  const std::vector<Acts::SourceLink> links{
      Acts::SourceLink{1}, Acts::SourceLink{2}, Acts::SourceLink{3}};
  for (bool extra : {false, true}) {
    for (bool outlier : {false, true}) {
      for (std::size_t count : {0u, 1u, 2u, 3u}) {
        Acts::VectorMultiTrajectory trajectory, reference;
        trajectory.addColumn<int>("marker");
        reference.addColumn<int>("marker");
        Creator creator;
        Calibrator calibrator{extra};
        Selector selector{count, outlier};
        creator.calibrator.connect<&Calibrator::calibrate>(&calibrator);
        creator.measurementSelector.connect<&Selector::select>(&selector);
        Creator::candidate_container_t candidates;
        auto result = creator.createSourceLinkTrackStates(
            gctx, cctx, *surface, bound, links.begin(), links.end(),
            Acts::kTrackIndexInvalid, candidates, trajectory,
            Acts::getDummyLogger());
        BOOST_REQUIRE(result.ok());
        auto copied = creator.processSelectedTrackStates(
            candidates.begin(), candidates.begin() + count, reference, outlier,
            Acts::getDummyLogger());
        BOOST_REQUIRE(copied.ok());
        BOOST_REQUIRE_EQUAL(result->size(), copied->size());
        for (std::size_t i = 0; i < count; ++i) {
          const auto a = trajectory.getTrackState(result->at(i));
          const auto b = reference.getTrackState(copied->at(i));
          BOOST_CHECK_EQUAL(a.previous(), b.previous());
          BOOST_CHECK_EQUAL(a.getUncalibratedSourceLink().get<int>(),
                            b.getUncalibratedSourceLink().get<int>());
          BOOST_CHECK_EQUAL(a.component<int>("marker"),
                            b.component<int>("marker"));
          BOOST_CHECK_EQUAL(a.pathLength(), b.pathLength());
          BOOST_CHECK_EQUAL(a.typeFlags().raw(), b.typeFlags().raw());
          BOOST_CHECK(!a.hasFiltered());
          for (std::size_t row = 0; row < Acts::eBoundSize; ++row) {
            BOOST_CHECK_EQUAL(a.predicted()[row], b.predicted()[row]);
            for (std::size_t col = 0; col < Acts::eBoundSize; ++col) {
              BOOST_CHECK_EQUAL(a.predictedCovariance()(row, col),
                                b.predictedCovariance()(row, col));
              BOOST_CHECK_EQUAL(a.jacobian()(row, col), b.jacobian()(row, col));
            }
          }
          for (int row = 0; row < 2; ++row) {
            BOOST_CHECK_EQUAL(a.calibrated<2>()[row], b.calibrated<2>()[row]);
            for (int col = 0; col < 2; ++col) {
              BOOST_CHECK_EQUAL(a.calibratedCovariance<2>()(row, col),
                                b.calibratedCovariance<2>()(row, col));
            }
          }
        }
      }
    }
  }
}
