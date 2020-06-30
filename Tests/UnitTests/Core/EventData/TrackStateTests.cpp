// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/optional/optional_io.hpp>
#include <boost/test/unit_test.hpp>

#include <random>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/EventData/TrackStateSorters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

using SourceLink = MinimalSourceLink;

template <ParID_t... params>
using MeasurementType = Measurement<SourceLink, params...>;
using FittableMeasurement = FittableMeasurement<SourceLink>;
using BoundTrackState = TrackState<SourceLink, BoundParameters>;
///
/// @brief Unit test for creation of Measurement object
///
BOOST_AUTO_TEST_CASE(track_state_initialization) {
  std::default_random_engine generator(42);

  // The plane surface
  auto plane = Surface::makeShared<PlaneSurface>(Vector3D{0., 0., 0.},
                                                 Vector3D{1., 0., 0.});

  // Construct the 1D measurement
  ActsSymMatrixD<1> cov1D;
  cov1D << 0.04;

  FittableMeasurement m1D(
      MeasurementType<ParDef::eLOC_0>(plane, {}, std::move(cov1D), 0.02));

  // Construct the 2D measurement
  SymMatrix2D cov2D;
  cov2D << 0.04, 0., 0.09, 0.;

  FittableMeasurement m2D(MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1>(
      plane, {}, std::move(cov2D), 0.02, 0.03));

  // The 1D track state from the measurement
  BoundTrackState mts1D(SourceLink{&m1D});

  auto calibrate = [](auto& trackstate) {
    // double ref: optional, MinimalSourceLink
    trackstate.measurement.calibrated = **trackstate.measurement.uncalibrated;
  };

  // "calibrate" the measurement
  calibrate(mts1D);

  BOOST_CHECK_EQUAL(*mts1D.size(), 1u);

  // Test the copy construtor
  BoundTrackState mts1DCopy(mts1D);

  // Test the copy move constructor
  BoundTrackState mts1DCopyMoved(std::move(mts1DCopy));

  // Test the copy assignment operator
  BoundTrackState mts1DCopyAssigned = mts1DCopyMoved;

  // Test the comovepy assignment operator
  BoundTrackState mts1DMoveAssigned = std::move(mts1DCopyAssigned);

  // Swap the measurements
  std::swap(mts1DMoveAssigned, mts1D);

  // The 2D track state from the measurement
  BoundTrackState mts2D(SourceLink{&m2D});
  calibrate(mts2D);

  BOOST_CHECK_EQUAL(*mts2D.size(), 2u);

  // Construct the parameter
  std::array<double, 6> pars_array = {
      {-0.1234, 9.8765, 0.45, 0.888, 0.001, 0.}};
  BoundVector pars;
  pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
      pars_array[4], pars_array[5];

  // constructor from parameter vector: predicted filtered, smoothed
  BoundParameters ataPlane(tgContext, std::nullopt, pars, plane);

  // The parameter track state from the parameters
  BoundTrackState pts(std::move(ataPlane));

  // Test the copy constructor for a parameter state
  BoundTrackState ptsCopy(pts);

  // Test the copy move constructor for a parameter state
  BoundTrackState ptsCopyMove(std::move(ptsCopy));

  // Test the copy assignment for a parameter state
  BoundTrackState ptsCopyAssigned = ptsCopyMove;

  // Test the move assignment for a parameter state
  BoundTrackState ptsMoveAssigned = std::move(ptsCopyAssigned);

  std::vector<BoundTrackState> trackStates = {std::move(mts1DMoveAssigned),
                                              std::move(mts2D), std::move(pts)};

  BOOST_CHECK_EQUAL(trackStates.size(), 3u);

  // Test is we can shuffle the track states
  // Test to extract the surface of these guys
  for (auto& ts : trackStates) {
    const Surface* sf = &ts.referenceSurface();
    BOOST_CHECK_EQUAL(sf, plane.get());
  }

  // Create predicted, filtered and smoothed parameters
  BoundParameters ataPlaneUpdt(tgContext, std::nullopt, pars, plane);
  BoundParameters ataPlanePred(tgContext, std::nullopt, pars, plane);
  BoundParameters ataPlaneSmth(tgContext, std::nullopt, pars, plane);

  // Get the predicted parameters back from the trackState
  auto& ptsfList = trackStates[2];
  auto& ataPlanefListPred = ptsfList.parameter.predicted;
  BOOST_CHECK(ataPlanefListPred);

  // Check that the other parameters are empty
  auto& ataPlanefListUpdt = ptsfList.parameter.filtered;
  BOOST_CHECK(!ataPlanefListUpdt);

  auto& ataPlanefListSmthd = ptsfList.parameter.smoothed;
  BOOST_CHECK(!ataPlanefListSmthd);

  // Get the track States from the list
  auto& m2DfList = trackStates[1];

  m2DfList.parameter.filtered = std::move(ataPlaneUpdt);
  auto& ataPlanefListUpdtM2D = m2DfList.parameter.filtered;
  BOOST_CHECK(ataPlanefListUpdtM2D);

  m2DfList.parameter.predicted = std::move(ataPlanePred);
  auto& ataPlanefListPred2D = m2DfList.parameter.predicted;
  BOOST_CHECK(ataPlanefListPred2D);

  // Test the sorting helper
  BoundParameters ataPlaneAt1(tgContext, std::nullopt, pars, plane);
  BoundTrackState ataPlaneState1(std::move(ataPlaneAt1));
  ataPlaneState1.parameter.pathLength = 1.;

  BoundParameters ataPlaneAt2(tgContext, std::nullopt, pars, plane);
  BoundTrackState ataPlaneState2(std::move(ataPlaneAt2));
  ataPlaneState2.parameter.pathLength = 2.;

  std::vector<BoundTrackState> unorderedStates = {std::move(ataPlaneState2),
                                                  std::move(ataPlaneState1)};

  // Sort the variant track state
  TrackStatePathLengthSorter plSorter;
  std::sort(unorderedStates.begin(), unorderedStates.end(), plSorter);

  auto firstOrdered = unorderedStates[0];
  BOOST_CHECK_EQUAL(firstOrdered.parameter.pathLength, 1.);

  auto secondOrdered = unorderedStates[1];
  BOOST_CHECK_EQUAL(secondOrdered.parameter.pathLength, 2.);

  auto& pState = firstOrdered.parameter;

  BOOST_CHECK_EQUAL(pState.pathLength, 1.);

  std::shuffle(unorderedStates.begin(), unorderedStates.end(), generator);

  // Copy the TrackStates into a new vector
  std::vector<BoundTrackState> copiedStates = {unorderedStates[0],
                                               unorderedStates[1]};

  // Shuffle
  std::shuffle(copiedStates.begin(), copiedStates.end(), generator);
  // And sort again
  std::sort(copiedStates.begin(), copiedStates.end(), plSorter);
}

}  // namespace Test
}  // namespace Acts
