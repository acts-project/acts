// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE TrackState Tests
#include <boost/test/included/unit_test.hpp>

#include <random>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/EventData/detail/surface_getter.hpp"
#include "Acts/EventData/detail/trackstate_manipulation.hpp"
#include "Acts/EventData/detail/trackstate_sorters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace Test {

  using Jacobian   = ActsMatrixD<5, 5>;
  using Identifier = unsigned long int;

  template <ParID_t... params>
  using MeasurementType = Measurement<Identifier, params...>;
  template <ParID_t... params>
  using MeasuredTrackState
      = MeasuredTrackState<Identifier, BoundParameters, Jacobian, params...>;
  using ParametricTrackState
      = ParametricTrackState<Identifier, BoundParameters, Jacobian>;
  using VariantState = VariantTrackState<Identifier, BoundParameters, Jacobian>;
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
    MeasuredTrackState<ParDef::eLOC_0> mts1D(std::move(m1D));

    // Test the copy construtor
    MeasuredTrackState<ParDef::eLOC_0> mts1DCopy(mts1D);

    // Test the copy move constructor
    MeasuredTrackState<ParDef::eLOC_0> mts1DCopyMoved(std::move(mts1DCopy));

    // Test the copy assignment operator
    MeasuredTrackState<ParDef::eLOC_0> mts1DCopyAssigned = mts1DCopyMoved;

    // Test the comovepy assignment operator
    MeasuredTrackState<ParDef::eLOC_0> mts1DMoveAssigned
        = std::move(mts1DCopyAssigned);

    // Swap the measurements
    std::swap(mts1DMoveAssigned, mts1D);

    VariantState m1DVariant(mts1D);

    VariantState m1DVariantMoveAssigned(mts1DMoveAssigned);

    std::swap(m1DVariant, m1DVariantMoveAssigned);

    // The 2D track state from the measurement
    MeasuredTrackState<ParDef::eLOC_0, ParDef::eLOC_1> mts2D(std::move(m2D));

    // Construct the parameter
    std::array<double, 5> pars_array = {{-0.1234, 9.8765, 0.45, 0.888, 0.001}};
    TrackParametersBase::ParVector_t pars;
    pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
        pars_array[4];

    // Make a variant copy of the 2D TrackState
    VariantState m2DVariant(mts2D);

    // constructor from parameter vector: predicted filtered, smoothed
    BoundParameters ataPlane(nullptr, pars, plane);

    // The parametric track state from the parameters
    ParametricTrackState pts(std::move(ataPlane));

    // Test the copy constructor for a parametric state
    ParametricTrackState ptsCopy(pts);

    // Test the copy move constructor for a parametric state
    ParametricTrackState ptsCopyMove(ptsCopy);

    // Test the copy assignment for a parametric state
    ParametricTrackState ptsCopyAssigned = ptsCopyMove;

    // Test the move assignment for a parametric state
    ParametricTrackState ptsMoveAssigned = std::move(ptsCopyAssigned);

    std::vector<VariantState> trackStates
        = {std::move(m1D), std::move(m2D), std::move(pts)};

    BOOST_CHECK(trackStates.size() == 3);

    std::vector<VariantState> unorderedStates{
        trackStates[2], trackStates[1], trackStates[0], trackStates[3]};

    // Test is we can shuffle the track states
    // Test to extract the surface of these guys
    for (auto& ts : trackStates) {
      const Surface* sf = &(detail::getSurface(ts));
      bool           surfaceExists(sf != nullptr);
      BOOST_TEST(surfaceExists);
    }
    // Create predicted, filtered and smoothed parameters
    BoundParameters ataPlaneUpdt(nullptr, pars, plane);
    BoundParameters ataPlanePred(nullptr, pars, plane);
    BoundParameters ataPlaneSmth(nullptr, pars, plane);

    // Get the predicted parameters back from the trackState
    auto& ptsfList          = trackStates[2];
    auto  ataPlanefListPred = detail::getParamaters<BoundParameters>(
        ptsfList, ParametricType::predicted);
    BOOST_TEST(ataPlanefListPred);

    // Check that the other parameters are empty
    auto ataPlanefListUpdt = detail::getParamaters<BoundParameters>(
        ptsfList, ParametricType::filtered);
    BOOST_TEST(!ataPlanefListUpdt);

    auto ataPlanefListSmthd = detail::getParamaters<BoundParameters>(
        ptsfList, ParametricType::smoothed);
    BOOST_TEST(!ataPlanefListSmthd);

    // Get the track States from the list
    auto& m2DfList = trackStates[1];

    detail::setParameters<BoundParameters>(
        m2DfList, std::move(ataPlaneUpdt), ParametricType::filtered);
    auto ataPlanefListUpdtM2D = detail::getParamaters<BoundParameters>(
        m2DfList, ParametricType::filtered);
    BOOST_TEST(ataPlanefListUpdtM2D);

    detail::setParameters<BoundParameters>(
        m2DfList, std::move(ataPlanePred), ParametricType::predicted);
    auto ataPlanefListPred2D = detail::getParamaters<BoundParameters>(
        m2DfList, ParametricType::predicted);
    BOOST_TEST(ataPlanefListPred2D);

    // Test the sorting helper
    ParametricTrackState surfacePlaneAt2(&plane);
    surfacePlaneAt2.parametric.pathLength = 2.;

    BoundParameters      ataPlaneAt1(nullptr, pars, plane);
    ParametricTrackState ataPlaneState1(std::move(ataPlaneAt1));
    ataPlaneState1.parametric.pathLength = 1.;

    std::vector<VariantState> undorderedStates
        = {std::move(surfacePlaneAt2), std::move(ataPlaneState1)};

    // Sort the variant track state
    detail::path_length_sorter plSorter;
    std::sort(undorderedStates.begin(), undorderedStates.end(), plSorter);

    auto   firstOrdered   = undorderedStates[0];
    double fistPathLength = detail::getPathLength(firstOrdered);
    BOOST_TEST(fistPathLength == 1.);

    auto   secondOrdered    = undorderedStates[1];
    double secondPathLength = detail::getPathLength(secondOrdered);
    BOOST_TEST(secondPathLength == 2.);

    using ParState = ParametricState<BoundParameters, Jacobian>;
    auto& pState   = detail::getParametricState<ParState>(firstOrdered);

    BOOST_TEST(pState.pathLength == 1.);

    std::shuffle(undorderedStates.begin(), undorderedStates.end(), generator);

    // Copy the TrackStates into a new vector
    std::vector<VariantState> copiedStates
        = {undorderedStates[0], undorderedStates[1]};

    // Shuffle
    std::shuffle(copiedStates.begin(), copiedStates.end(), generator);
    // And sort again
    std::sort(copiedStates.begin(), copiedStates.end(), plSorter);
  }

}  // namespace Test
}  // namespace Acts
