// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE TrackState Tests
#include <boost/optional/optional_io.hpp>
#include <boost/test/included/unit_test.hpp>

#include <random>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/EventData/TrackStateSorters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace Test {

  using Jacobian   = ActsMatrixD<5, 5>;
  using Identifier = unsigned long int;

  template <ParID_t... params>
  using MeasurementType = Measurement<Identifier, params...>;
  using BoundTrackState = TrackState<Identifier, BoundParameters>;
  ///
  /// @brief Unit test for creation of Measurement object
  ///
  BOOST_AUTO_TEST_CASE(track_state_initialization)
  {

    std::default_random_engine generator(42);

    // The plane surface
    auto plane = Surface::makeShared<PlaneSurface>(Vector3D{0., 0., 0.},
                                                   Vector3D{1., 0., 0.});

    // Construct the 1D measurement
    ActsSymMatrixD<1> cov1D;
    cov1D << 0.04;

    MeasurementType<ParDef::eLOC_0> m1D(plane, 0, std::move(cov1D), 0.02);

    // Construct the 2D measurement
    ActsSymMatrixD<2> cov2D;
    cov2D << 0.04, 0., 0.09, 0.;

    MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1> m2D(
        plane, 0, std::move(cov2D), 0.02, 0.03);

    // The 1D track state from the measurement
    BoundTrackState mts1D(std::move(m1D));

    BOOST_CHECK_EQUAL(*mts1D.size(), 1);

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
    BoundTrackState mts2D(std::move(m2D));

    BOOST_CHECK_EQUAL(*mts2D.size(), 2);

    // Construct the parameter
    std::array<double, 5> pars_array = {{-0.1234, 9.8765, 0.45, 0.888, 0.001}};
    TrackParametersBase::ParVector_t pars;
    pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
        pars_array[4];

    // constructor from parameter vector: predicted filtered, smoothed
    BoundParameters ataPlane(nullptr, pars, plane);

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

    std::vector<BoundTrackState> trackStates
        = {std::move(mts1DMoveAssigned), std::move(mts2D), std::move(pts)};

    BOOST_CHECK_EQUAL(trackStates.size(), 3);

    // Test is we can shuffle the track states
    // Test to extract the surface of these guys
    for (auto& ts : trackStates) {
      const Surface* sf = &ts.referenceSurface();
      BOOST_CHECK_EQUAL(sf, plane.get());
    }

    // Create predicted, filtered and smoothed parameters
    BoundParameters ataPlaneUpdt(nullptr, pars, plane);
    BoundParameters ataPlanePred(nullptr, pars, plane);
    BoundParameters ataPlaneSmth(nullptr, pars, plane);

    // Get the predicted parameters back from the trackState
    auto& ptsfList          = trackStates[2];
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
    auto& ataPlanefListUpdtM2D  = m2DfList.parameter.filtered;
    BOOST_CHECK(ataPlanefListUpdtM2D);

    m2DfList.parameter.predicted = std::move(ataPlanePred);
    auto& ataPlanefListPred2D    = m2DfList.parameter.predicted;
    BOOST_CHECK(ataPlanefListPred2D);

    // Test the sorting helper
    BoundParameters ataPlaneAt1(nullptr, pars, plane);
    BoundTrackState ataPlaneState1(std::move(ataPlaneAt1));
    ataPlaneState1.parameter.pathLength = 1.;

    BoundParameters ataPlaneAt2(nullptr, pars, plane);
    BoundTrackState ataPlaneState2(std::move(ataPlaneAt2));
    ataPlaneState2.parameter.pathLength = 2.;

    std::vector<BoundTrackState> unorderedStates
        = {std::move(ataPlaneState2), std::move(ataPlaneState1)};

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
    std::vector<BoundTrackState> copiedStates
        = {unorderedStates[0], unorderedStates[1]};

    // Shuffle
    std::shuffle(copiedStates.begin(), copiedStates.end(), generator);
    // And sort again
    std::sort(copiedStates.begin(), copiedStates.end(), plSorter);
  }

}  // namespace Test
}  // namespace Acts
