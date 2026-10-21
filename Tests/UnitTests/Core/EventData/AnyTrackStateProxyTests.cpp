// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackStateProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/HashedString.hpp"

namespace {

using namespace Acts;
using namespace Acts::HashedStringLiteral;

struct TestTrackStateFixture {
  using Trajectory = VectorMultiTrajectory;
  using TrackContainerBackend = VectorTrackContainer;
  using Container =
      TrackContainer<TrackContainerBackend, Trajectory, detail::RefHolder>;

  TestTrackStateFixture() : container(trackContainer, trajectory) {}

  Container container;
  TrackContainerBackend trackContainer;
  Trajectory trajectory;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataAnyTrackState)

BOOST_FIXTURE_TEST_CASE(WrapTrackStateProxy, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  auto proxy = container.trackStateContainer().getTrackState(state.index());
  AnyMutableTrackStateProxy anyState(proxy);

  BOOST_CHECK_EQUAL(anyState.index(), proxy.index());
}

BOOST_FIXTURE_TEST_CASE(ConstructFromReadOnlyTrackStateProxy,
                        TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();
  state.predicted() = Vector<eBoundSize>::Ones();

  TrackContainer constContainer{ConstVectorTrackContainer{trackContainer},
                                ConstVectorMultiTrajectory{trajectory}};

  auto constState =
      constContainer.trackStateContainer().getTrackState(state.index());
  AnyConstTrackStateProxy anyState(constState);

  BOOST_CHECK_EQUAL(anyState.index(), state.index());
  BOOST_CHECK_CLOSE(anyState.predicted()[eBoundLoc0], 1.0, 1e-6);
}

BOOST_FIXTURE_TEST_CASE(AccessFiltered, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  state.predicted() = Vector<eBoundSize>::Ones();
  state.filtered() = 2. * Vector<eBoundSize>::Ones();

  AnyMutableTrackStateProxy anyState(state);

  BOOST_CHECK_EQUAL(anyState.predicted()[eBoundLoc0], 1.);
  BOOST_CHECK_EQUAL(anyState.filtered()[eBoundLoc0], 2.);

  auto filtered = anyState.filtered();
  filtered[eBoundLoc0] = 3.;
  BOOST_CHECK_EQUAL(state.filtered()[eBoundLoc0], 3.);
}

BOOST_FIXTURE_TEST_CASE(AccessCalibratedFixedSize, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  state.allocateCalibrated(2);
  auto calibrated = state.template calibrated<2>();
  calibrated[0] = 1.5;
  calibrated[1] = 2.5;
  auto calibratedCov = state.template calibratedCovariance<2>();
  calibratedCov.setZero();
  calibratedCov(0, 0) = 0.1;
  calibratedCov(1, 1) = 0.2;

  AnyMutableTrackStateProxy anyState(state);

  auto view = anyState.calibrated<2>();
  BOOST_CHECK_CLOSE(view[0], 1.5, 1e-6);
  BOOST_CHECK_CLOSE(view[1], 2.5, 1e-6);

  view[0] = 4.5;
  BOOST_CHECK_CLOSE(state.calibrated<2>()[0], 4.5, 1e-6);

  auto covView = anyState.calibratedCovariance<2>();
  BOOST_CHECK_CLOSE(covView(0, 0), 0.1, 1e-6);
  BOOST_CHECK_CLOSE(covView(1, 1), 0.2, 1e-6);
  covView(0, 0) = 0.5;
  BOOST_CHECK_CLOSE(state.calibratedCovariance<2>()(0, 0), 0.5, 1e-6);

  AnyConstTrackStateProxy constState(state);
  auto constView = constState.calibrated<2>();
  BOOST_CHECK_CLOSE(constView[0], 4.5, 1e-6);
  auto constCov = constState.calibratedCovariance<2>();
  BOOST_CHECK_CLOSE(constCov(0, 0), 0.5, 1e-6);
}

BOOST_FIXTURE_TEST_CASE(AccessEffectiveCalibratedDynamic,
                        TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  constexpr std::size_t measdim = 3u;
  state.allocateCalibrated(measdim);
  {
    auto dyn = state.effectiveCalibrated();
    for (std::size_t i = 0; i < measdim; ++i) {
      dyn[i] = 1. + static_cast<double>(i);
    }
  }
  {
    auto dynCov = state.effectiveCalibratedCovariance();
    dynCov.setZero();
    for (std::size_t i = 0; i < measdim; ++i) {
      dynCov(i, i) = 0.1 * static_cast<double>(i + 1);
    }
  }

  AnyMutableTrackStateProxy anyState(state);
  BOOST_CHECK_EQUAL(anyState.calibratedSize(), measdim);
  auto eff = anyState.effectiveCalibrated();
  BOOST_CHECK_EQUAL(eff.size(), measdim);
  BOOST_CHECK_CLOSE(eff[0], 1., 1e-6);
  BOOST_CHECK_CLOSE(eff[2], 3., 1e-6);
  eff[1] = 5.5;
  BOOST_CHECK_CLOSE(state.effectiveCalibrated()[1], 5.5, 1e-6);

  auto effCov = anyState.effectiveCalibratedCovariance();
  BOOST_CHECK_EQUAL(effCov.rows(), measdim);
  BOOST_CHECK_EQUAL(effCov.cols(), measdim);
  BOOST_CHECK_CLOSE(effCov(0, 0), 0.1, 1e-6);
  BOOST_CHECK_CLOSE(effCov(2, 2), 0.3, 1e-6);
  effCov(1, 1) = 1.7;
  BOOST_CHECK_CLOSE(state.effectiveCalibratedCovariance()(1, 1), 1.7, 1e-6);

  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK_EQUAL(constState.calibratedSize(), measdim);
  auto constEff = constState.effectiveCalibrated();
  BOOST_CHECK_EQUAL(constEff.size(), measdim);
  BOOST_CHECK_CLOSE(constEff[1], 5.5, 1e-6);
  auto constEffCov = constState.effectiveCalibratedCovariance();
  BOOST_CHECK_EQUAL(constEffCov.rows(), measdim);
  BOOST_CHECK_CLOSE(constEffCov(1, 1), 1.7, 1e-6);
}

BOOST_FIXTURE_TEST_CASE(ProxyAccessorWithAnyTrackState, TestTrackStateFixture) {
  container.trackStateContainer().addColumn<float>("customFloat");

  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  state.template component<float>("customFloat"_hash) = 0.25f;

  ProxyAccessor<float> mutableAccessor("customFloat");
  ConstProxyAccessor<float> constAccessor("customFloat");

  AnyMutableTrackStateProxy anyState(state);
  BOOST_CHECK_CLOSE(mutableAccessor(anyState), 0.25f, 1e-6);
  mutableAccessor(anyState) = 0.75f;
  BOOST_CHECK_CLOSE(state.template component<float>("customFloat"_hash), 0.75f,
                    1e-6);

  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK_CLOSE(constAccessor(constState), 0.75f, 1e-6);
  BOOST_CHECK(constAccessor.hasColumn(constState));
}

BOOST_FIXTURE_TEST_CASE(PreviousAndHasPrevious, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state1 = container.trackStateContainer().makeTrackState();
  auto state2 = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state2.index();

  // Link state2 -> state1
  state2.previous() = state1.index();

  AnyMutableTrackStateProxy anyState1(state1);
  AnyMutableTrackStateProxy anyState2(state2);

  // state1 has no previous
  BOOST_CHECK(!anyState1.hasPrevious());

  // state2 has previous
  BOOST_CHECK(anyState2.hasPrevious());
  BOOST_CHECK_EQUAL(anyState2.previous(), state1.index());

  // Test with const proxy
  AnyConstTrackStateProxy constState2(state2);
  BOOST_CHECK(constState2.hasPrevious());
  BOOST_CHECK_EQUAL(constState2.previous(), state1.index());
}

BOOST_FIXTURE_TEST_CASE(GetMask, TestTrackStateFixture) {
  // Test with no components
  {
    auto state = container.trackStateContainer().makeTrackState(
        TrackStatePropMask::None);
    AnyMutableTrackStateProxy anyState(state);
    auto mask = anyState.getMask();
    BOOST_CHECK_EQUAL(mask, TrackStatePropMask::None);
  }

  // Test with only predicted
  {
    auto state = container.trackStateContainer().makeTrackState(
        TrackStatePropMask::Predicted);
    AnyMutableTrackStateProxy anyState(state);
    auto mask = anyState.getMask();
    BOOST_CHECK(ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted));
    BOOST_CHECK(!ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered));
    BOOST_CHECK(!ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed));
  }

  // Test with predicted and filtered
  {
    auto state = container.trackStateContainer().makeTrackState(
        TrackStatePropMask::Predicted | TrackStatePropMask::Filtered);
    AnyMutableTrackStateProxy anyState(state);
    auto mask = anyState.getMask();
    BOOST_CHECK(ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted));
    BOOST_CHECK(ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered));
    BOOST_CHECK(!ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed));
  }

  // Test with all components
  {
    auto state =
        container.trackStateContainer().makeTrackState(TrackStatePropMask::All);
    state.allocateCalibrated(2);
    AnyMutableTrackStateProxy anyState(state);
    auto mask = anyState.getMask();
    BOOST_CHECK(ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted));
    BOOST_CHECK(ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered));
    BOOST_CHECK(ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed));
    BOOST_CHECK(ACTS_CHECK_BIT(mask, TrackStatePropMask::Jacobian));
    BOOST_CHECK(ACTS_CHECK_BIT(mask, TrackStatePropMask::Calibrated));
  }
}

BOOST_FIXTURE_TEST_CASE(AccessSmoothed, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  // Set smoothed parameters and covariance
  state.smoothed() = Vector<eBoundSize>::Constant(5.);
  state.smoothedCovariance() = Matrix<eBoundSize, eBoundSize>::Identity();

  AnyMutableTrackStateProxy anyState(state);

  // Check hasSmoothed
  BOOST_CHECK(anyState.hasSmoothed());

  // Access smoothed parameters
  auto smoothed = anyState.smoothed();
  BOOST_CHECK_CLOSE(smoothed[eBoundLoc0], 5., 1e-6);

  // Modify smoothed parameters
  smoothed[eBoundLoc0] = 7.;
  BOOST_CHECK_CLOSE(state.smoothed()[eBoundLoc0], 7., 1e-6);

  // Access smoothed covariance
  auto smoothedCov = anyState.smoothedCovariance();
  BOOST_CHECK_CLOSE(smoothedCov(0, 0), 1., 1e-6);
  smoothedCov(0, 0) = 2.;
  BOOST_CHECK_CLOSE(state.smoothedCovariance()(0, 0), 2., 1e-6);

  // Test const proxy
  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK(constState.hasSmoothed());
  auto constSmoothed = constState.smoothed();
  BOOST_CHECK_CLOSE(constSmoothed[eBoundLoc0], 7., 1e-6);
}

BOOST_FIXTURE_TEST_CASE(AccessJacobian, TestTrackStateFixture) {
  auto track = container.makeTrack();
  // Create state with jacobian allocated
  auto state = container.trackStateContainer().makeTrackState(
      TrackStatePropMask::Jacobian);
  track.tipIndex() = state.index();

  AnyMutableTrackStateProxy anyState(state);
  BOOST_CHECK(anyState.hasJacobian());

  // Set jacobian values
  state.jacobian() = Matrix<eBoundSize, eBoundSize>::Identity();
  state.jacobian()(0, 1) = 0.5;

  // Access jacobian (mutable)
  auto jac = anyState.jacobian();
  BOOST_CHECK_CLOSE(jac(0, 0), 1., 1e-6);
  BOOST_CHECK_CLOSE(jac(0, 1), 0.5, 1e-6);

  // Modify jacobian
  jac(0, 1) = 0.75;
  BOOST_CHECK_CLOSE(state.jacobian()(0, 1), 0.75, 1e-6);

  // Test const proxy
  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK(constState.hasJacobian());
  auto constJac = constState.jacobian();
  BOOST_CHECK_CLOSE(constJac(0, 1), 0.75, 1e-6);
}

BOOST_FIXTURE_TEST_CASE(AccessParametersAndCovariance, TestTrackStateFixture) {
  auto track = container.makeTrack();
  // Allocate predicted parameters
  auto state = container.trackStateContainer().makeTrackState(
      TrackStatePropMask::Predicted);
  track.tipIndex() = state.index();

  // Set predicted parameters and covariance
  state.predicted() = Vector<eBoundSize>::Constant(3.);
  state.predictedCovariance() = Matrix<eBoundSize, eBoundSize>::Identity();
  state.predictedCovariance()(1, 1) = 2.;

  AnyMutableTrackStateProxy anyState(state);

  // Access through parameters() and covariance() directly
  auto params = anyState.parameters();
  BOOST_CHECK_CLOSE(params[eBoundLoc0], 3., 1e-6);

  auto cov = anyState.covariance();
  BOOST_CHECK_CLOSE(cov(0, 0), 1., 1e-6);
  BOOST_CHECK_CLOSE(cov(1, 1), 2., 1e-6);

  // Access through predictedCovariance
  auto predCov = anyState.predictedCovariance();
  BOOST_CHECK_CLOSE(predCov(1, 1), 2., 1e-6);
  predCov(1, 1) = 3.;
  BOOST_CHECK_CLOSE(state.predictedCovariance()(1, 1), 3., 1e-6);

  // Test const proxy
  AnyConstTrackStateProxy constState(state);
  auto constParams = constState.parameters();
  BOOST_CHECK_CLOSE(constParams[eBoundLoc0], 3., 1e-6);
  auto constCov = constState.covariance();
  BOOST_CHECK_CLOSE(constCov(1, 1), 3., 1e-6);
}

BOOST_FIXTURE_TEST_CASE(ReferenceSurface, TestTrackStateFixture) {
  auto track = container.makeTrack();
  // Create with default mask which includes all components
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  // Create and set a reference surface
  auto initialSurface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(10, 10));
  state.setReferenceSurface(initialSurface);

  AnyMutableTrackStateProxy anyState(state);

  // Should now have a reference surface
  BOOST_CHECK(anyState.hasReferenceSurface());
  const Surface& surf1 = anyState.referenceSurface();
  BOOST_CHECK(&surf1 == initialSurface.get());

  // Create a new surface and set it
  auto newSurface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(10, 10));
  anyState.setReferenceSurface(newSurface);

  BOOST_CHECK(anyState.hasReferenceSurface());
  const Surface& surf2 = anyState.referenceSurface();
  BOOST_CHECK(&surf2 == newSurface.get());

  // Test const proxy
  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK(constState.hasReferenceSurface());
  const Surface& surf3 = constState.referenceSurface();
  BOOST_CHECK(&surf3 == newSurface.get());
}

BOOST_FIXTURE_TEST_CASE(SourceLink, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  AnyMutableTrackStateProxy anyState(state);

  // Initially should not have source link
  BOOST_CHECK(!anyState.hasUncalibratedSourceLink());

  // Set source link
  state.setUncalibratedSourceLink(Acts::SourceLink{42});

  // Now should have source link
  BOOST_CHECK(anyState.hasUncalibratedSourceLink());

  // Get source link
  auto retrievedLink = anyState.getUncalibratedSourceLink();
  BOOST_CHECK_EQUAL(retrievedLink.get<int>(), 42);

  // Set new source link through any proxy
  Acts::SourceLink newLink{99};
  anyState.setUncalibratedSourceLink(std::move(newLink));

  // Should still have source link
  BOOST_CHECK(anyState.hasUncalibratedSourceLink());

  auto retrievedLink2 = anyState.getUncalibratedSourceLink();
  BOOST_CHECK_EQUAL(retrievedLink2.get<int>(), 99);

  // Test const proxy
  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK(constState.hasUncalibratedSourceLink());
  auto constLink = constState.getUncalibratedSourceLink();
  BOOST_CHECK_EQUAL(constLink.get<int>(), 99);
}

BOOST_FIXTURE_TEST_CASE(ProjectorSubspaceIndices, TestTrackStateFixture) {
  auto track = container.makeTrack();
  // Create with all components (projector is allocated with calibrated)
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  // Allocate calibrated which may also allocate projector
  constexpr std::size_t measdim = 2;
  state.allocateCalibrated(measdim);

  AnyMutableTrackStateProxy anyState(state);
  // Projector is now allocated
  BOOST_CHECK(anyState.hasProjector());

  // Set projector subspace indices
  BoundSubspaceIndices indices{};
  indices[0] = eBoundLoc0;
  indices[1] = eBoundLoc1;
  state.setProjectorSubspaceIndices(
      std::span<std::uint8_t>(indices.begin(), measdim));

  BOOST_CHECK(anyState.hasProjector());

  // Get projector subspace indices
  auto retrieved = anyState.projectorSubspaceIndices();
  BOOST_CHECK_EQUAL(retrieved[0], eBoundLoc0);
  BOOST_CHECK_EQUAL(retrieved[1], eBoundLoc1);

  // Get fixed-size variant
  auto fixedIndices = anyState.projectorSubspaceIndices<measdim>();
  BOOST_CHECK_EQUAL(fixedIndices[0], eBoundLoc0);
  BOOST_CHECK_EQUAL(fixedIndices[1], eBoundLoc1);

  // Test subspace helpers
  auto variableHelper = anyState.projectorSubspaceHelper();
  BOOST_CHECK_EQUAL(variableHelper.size(), measdim);

  auto fixedHelper = anyState.projectorSubspaceHelper<measdim>();
  BOOST_CHECK_EQUAL(fixedHelper.size(), measdim);

  // Test setting through any proxy
  BoundSubspaceIndices newIndices{};
  newIndices[0] = eBoundPhi;
  newIndices[1] = eBoundTheta;
  anyState.setProjectorSubspaceIndices(
      std::span<std::uint8_t>(newIndices.begin(), measdim));

  auto retrieved2 = anyState.projectorSubspaceIndices();
  BOOST_CHECK_EQUAL(retrieved2[0], eBoundPhi);
  BOOST_CHECK_EQUAL(retrieved2[1], eBoundTheta);

  // Test const proxy
  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK(constState.hasProjector());
  auto constIndices = constState.projectorSubspaceIndices();
  BOOST_CHECK_EQUAL(constIndices[0], eBoundPhi);
}

BOOST_FIXTURE_TEST_CASE(Chi2AndPathLength, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  // Set chi2 and pathLength
  state.chi2() = 2.5f;
  state.pathLength() = 123.456;

  AnyMutableTrackStateProxy anyState(state);

  // Access chi2
  BOOST_CHECK_CLOSE(anyState.chi2(), 2.5f, 1e-6);

  // Modify chi2
  anyState.chi2() = 3.5f;
  BOOST_CHECK_CLOSE(state.chi2(), 3.5f, 1e-6);

  // Access pathLength
  BOOST_CHECK_CLOSE(anyState.pathLength(), 123.456, 1e-6);

  // Modify pathLength
  anyState.pathLength() = 234.567;
  BOOST_CHECK_CLOSE(state.pathLength(), 234.567, 1e-6);

  // Test const proxy
  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK_CLOSE(constState.chi2(), 3.5f, 1e-6);
  BOOST_CHECK_CLOSE(constState.pathLength(), 234.567, 1e-6);
}

BOOST_FIXTURE_TEST_CASE(TypeFlags, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  // Initially typeFlags should be accessible (always allocated)
  AnyMutableTrackStateProxy anyState(state);

  // Set type flags through the original proxy
  state.typeFlags().setIsOutlier();

  // Access typeFlags through AnyTrackStateProxy
  auto flags = anyState.typeFlags();
  BOOST_CHECK(flags.hasMeasurement());
  BOOST_CHECK(!flags.isMeasurement());
  BOOST_CHECK(flags.isOutlier());
  BOOST_CHECK(!flags.isHole());

  // Modify typeFlags through AnyTrackStateProxy
  anyState.typeFlags().setIsHole();

  // Access typeFlags through AnyTrackStateProxy
  flags = anyState.typeFlags();
  BOOST_CHECK(!flags.hasMeasurement());
  BOOST_CHECK(!flags.isMeasurement());
  BOOST_CHECK(!flags.isOutlier());
  BOOST_CHECK(flags.isHole());
}

BOOST_FIXTURE_TEST_CASE(ComponentAccessWithStringKeys, TestTrackStateFixture) {
  container.trackStateContainer().addColumn<int>("intColumn");
  container.trackStateContainer().addColumn<double>("doubleColumn");

  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  state.template component<int>("intColumn"_hash) = 42;
  state.template component<double>("doubleColumn"_hash) = 3.14;

  AnyMutableTrackStateProxy anyState(state);

  // Test has() with string_view
  BOOST_CHECK(anyState.has("intColumn"));
  BOOST_CHECK(anyState.has("doubleColumn"));
  BOOST_CHECK(!anyState.has("nonexistent"));

  // Test component() with string_view (const)
  const auto& intVal = anyState.component<int>("intColumn");
  BOOST_CHECK_EQUAL(intVal, 42);
  const auto& doubleVal = anyState.component<double>("doubleColumn");
  BOOST_CHECK_CLOSE(doubleVal, 3.14, 1e-6);

  // Test component() with string_view (mutable)
  anyState.component<int>("intColumn") = 99;
  BOOST_CHECK_EQUAL(state.template component<int>("intColumn"_hash), 99);

  // Test hasColumn() with HashedString and string_view
  BOOST_CHECK(anyState.hasColumn("intColumn"_hash));
  BOOST_CHECK(anyState.hasColumn("intColumn"));
  BOOST_CHECK(!anyState.hasColumn("nonexistent"));

  // Test const proxy
  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK(constState.has("intColumn"));
  const auto& constIntVal = constState.component<int>("intColumn");
  BOOST_CHECK_EQUAL(constIntVal, 99);
}

BOOST_FIXTURE_TEST_CASE(UnsetComponents, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  // Add various components
  state.predicted() = Vector<eBoundSize>::Ones();
  state.filtered() = 2. * Vector<eBoundSize>::Ones();
  state.smoothed() = 3. * Vector<eBoundSize>::Ones();
  state.jacobian() = Matrix<eBoundSize, eBoundSize>::Identity();
  state.allocateCalibrated(2);

  AnyMutableTrackStateProxy anyState(state);

  // Verify all components exist
  BOOST_CHECK(anyState.hasPredicted());
  BOOST_CHECK(anyState.hasFiltered());
  BOOST_CHECK(anyState.hasSmoothed());
  BOOST_CHECK(anyState.hasJacobian());
  BOOST_CHECK(anyState.hasCalibrated());

  state.unset(TrackStatePropMask::Smoothed);

  // Verify through AnyTrackStateProxy that the component was unset
  BOOST_CHECK(anyState.hasPredicted());
  BOOST_CHECK(anyState.hasFiltered());
  BOOST_CHECK(!anyState.hasSmoothed());
  BOOST_CHECK(anyState.hasJacobian());
  BOOST_CHECK(anyState.hasCalibrated());

  anyState.unset(TrackStatePropMask::Jacobian);
  BOOST_CHECK(state.hasPredicted());
  BOOST_CHECK(state.hasFiltered());
  BOOST_CHECK(!state.hasSmoothed());
  BOOST_CHECK(!state.hasJacobian());
  BOOST_CHECK(state.hasCalibrated());

  anyState.unset(TrackStatePropMask::Calibrated);
  BOOST_CHECK(state.hasPredicted());
  BOOST_CHECK(state.hasFiltered());
  BOOST_CHECK(!anyState.hasSmoothed());
  BOOST_CHECK(!anyState.hasJacobian());
  BOOST_CHECK(!anyState.hasCalibrated());
}

BOOST_FIXTURE_TEST_CASE(AllocateCalibratedWithEigen, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  AnyMutableTrackStateProxy anyState(state);

  // Allocate with Eigen objects
  Vector<2> val;
  val << 1.5, 2.5;
  SquareMatrix<2> cov;
  cov << 0.1, 0.01, 0.01, 0.2;

  anyState.allocateCalibrated(val, cov);

  // Verify allocation
  BOOST_CHECK(anyState.hasCalibrated());
  BOOST_CHECK_EQUAL(anyState.calibratedSize(), 2);

  // Verify values
  auto calibrated = anyState.calibrated<2>();
  BOOST_CHECK_CLOSE(calibrated[0], 1.5, 1e-6);
  BOOST_CHECK_CLOSE(calibrated[1], 2.5, 1e-6);

  // Verify covariance
  auto calibratedCov = anyState.calibratedCovariance<2>();
  BOOST_CHECK_CLOSE(calibratedCov(0, 0), 0.1, 1e-6);
  BOOST_CHECK_CLOSE(calibratedCov(0, 1), 0.01, 1e-6);
  BOOST_CHECK_CLOSE(calibratedCov(1, 0), 0.01, 1e-6);
  BOOST_CHECK_CLOSE(calibratedCov(1, 1), 0.2, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
