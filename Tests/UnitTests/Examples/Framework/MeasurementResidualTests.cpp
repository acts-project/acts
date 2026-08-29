// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Exercises `measurementResidual`, which compares track state parameters
// against the state's own calibrated measurement.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Validation/ParametersOnSurface.hpp"

#include <array>
#include <memory>
#include <vector>

using namespace Acts;
using namespace ActsExamples;

namespace {

/// Local coordinates of the measurement, chosen so that both components of a
/// 2D measurement differ from the parameters below.
constexpr double measLoc0 = 1.5;
constexpr double measLoc1 = -2.5;
/// Variances of the measurement.
constexpr double measVar0 = 0.04;
constexpr double measVar1 = 0.09;
/// Local coordinates of the reconstructed parameters.
constexpr double recoLoc0 = 1.7;
constexpr double recoLoc1 = -2.1;
/// Variances of the reconstructed parameters.
constexpr double recoVar0 = 0.01;
constexpr double recoVar1 = 0.02;

std::shared_ptr<PlaneSurface> makeSurface() {
  return Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(100., 100.));
}

/// Reconstructed parameters with a diagonal covariance on @p surface.
BoundTrackParameters makeParameters(
    const std::shared_ptr<const Surface>& surface) {
  BoundVector parameters = BoundVector::Zero();
  parameters[eBoundLoc0] = recoLoc0;
  parameters[eBoundLoc1] = recoLoc1;

  BoundMatrix covariance = BoundMatrix::Zero();
  covariance(eBoundLoc0, eBoundLoc0) = recoVar0;
  covariance(eBoundLoc1, eBoundLoc1) = recoVar1;

  return BoundTrackParameters(surface, parameters, covariance,
                              ParticleHypothesis::pion());
}

/// Holds a one-state track container so the proxies stay valid.
struct Fixture {
  VectorTrackContainer trackBackend;
  VectorMultiTrajectory trajectory;
  Acts::TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                       Acts::detail::RefHolder>
      container{trackBackend, trajectory};

  std::shared_ptr<PlaneSurface> surface = makeSurface();

  /// Add a state carrying a measurement of @p indices.
  ///
  /// @param indices the bound indices the measurement constrains
  /// @return the state as a const proxy
  ConstTrackStateProxy addState(std::span<const std::uint8_t> indices) {
    auto state = container.trackStateContainer().makeTrackState();
    state.setReferenceSurface(surface);
    state.setProjectorSubspaceIndices(indices);

    state.allocateCalibrated(indices.size());
    auto calibrated = state.effectiveCalibrated();
    auto calibratedCovariance = state.effectiveCalibratedCovariance();
    calibratedCovariance.setZero();

    const std::array<double, 2> values{measLoc0, measLoc1};
    const std::array<double, 2> variances{measVar0, measVar1};
    for (std::size_t i = 0; i < indices.size(); ++i) {
      calibrated[i] = values.at(i);
      calibratedCovariance(i, i) = variances.at(i);
    }

    constContainer = std::make_unique<ConstTrackContainer>(
        std::make_shared<ConstVectorTrackContainer>(trackBackend),
        std::make_shared<ConstVectorMultiTrajectory>(trajectory));
    return constContainer->trackStateContainer().getTrackState(state.index());
  }

  std::unique_ptr<ConstTrackContainer> constContainer;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(ValidationMeasurementResidual)

BOOST_AUTO_TEST_CASE(UseOwnMeasurement) {
  using enum TrackParameterType;

  // only the parameters that consumed the state's own measurement are biased
  // towards it
  BOOST_CHECK(!parametersUseOwnMeasurement(Predicted));
  BOOST_CHECK(parametersUseOwnMeasurement(Filtered));
  BOOST_CHECK(parametersUseOwnMeasurement(Smoothed));
  BOOST_CHECK(!parametersUseOwnMeasurement(Unbiased));
}

BOOST_AUTO_TEST_CASE(TwoDimensionalMeasurement) {
  Fixture fixture;
  const std::array<std::uint8_t, 2> indices{eBoundLoc0, eBoundLoc1};
  const ConstTrackStateProxy state = fixture.addState(indices);

  const auto residual = measurementResidual(
      state, makeParameters(fixture.surface), TrackParameterType::Predicted);
  BOOST_REQUIRE(residual.has_value());

  BOOST_CHECK_EQUAL(residual->subspace.size(), 2u);
  BOOST_CHECK(residual->subspace.contains(eBoundLoc0));
  BOOST_CHECK(residual->subspace.contains(eBoundLoc1));

  // the reference is subtracted from the reconstructed value
  BOOST_CHECK_CLOSE(residual->residual[eBoundLoc0], recoLoc0 - measLoc0, 1e-9);
  BOOST_CHECK_CLOSE(residual->residual[eBoundLoc1], recoLoc1 - measLoc1, 1e-9);

  // predicted parameters do not use the measurement, so the covariances add
  BOOST_CHECK_CLOSE(residual->covariance(eBoundLoc0, eBoundLoc0),
                    measVar0 + recoVar0, 1e-9);
  BOOST_CHECK_CLOSE(residual->covariance(eBoundLoc1, eBoundLoc1),
                    measVar1 + recoVar1, 1e-9);

  // everything outside the measured subspace stays zero
  BOOST_CHECK_EQUAL(residual->residual[eBoundQOverP], 0.);
  BOOST_CHECK_EQUAL(residual->covariance(eBoundTheta, eBoundTheta), 0.);
}

BOOST_AUTO_TEST_CASE(SmoothedParametersSubtractTheCovariance) {
  Fixture fixture;
  const std::array<std::uint8_t, 2> indices{eBoundLoc0, eBoundLoc1};
  const ConstTrackStateProxy state = fixture.addState(indices);

  const auto residual = measurementResidual(
      state, makeParameters(fixture.surface), TrackParameterType::Smoothed);
  BOOST_REQUIRE(residual.has_value());

  // smoothed parameters used the measurement and are correlated with it
  BOOST_CHECK_CLOSE(residual->covariance(eBoundLoc0, eBoundLoc0),
                    measVar0 - recoVar0, 1e-9);
  BOOST_CHECK_CLOSE(residual->covariance(eBoundLoc1, eBoundLoc1),
                    measVar1 - recoVar1, 1e-9);

  // the residual itself does not depend on the parameter type
  BOOST_CHECK_CLOSE(residual->residual[eBoundLoc0], recoLoc0 - measLoc0, 1e-9);
}

BOOST_AUTO_TEST_CASE(OneDimensionalMeasurement) {
  Fixture fixture;
  const std::array<std::uint8_t, 1> indices{eBoundLoc0};
  const ConstTrackStateProxy state = fixture.addState(indices);

  const auto residual = measurementResidual(
      state, makeParameters(fixture.surface), TrackParameterType::Predicted);
  BOOST_REQUIRE(residual.has_value());

  // a strip constrains loc0 only, so loc1 is not part of the subspace and
  // must not be filled by the caller
  BOOST_CHECK_EQUAL(residual->subspace.size(), 1u);
  BOOST_CHECK(residual->subspace.contains(eBoundLoc0));
  BOOST_CHECK(!residual->subspace.contains(eBoundLoc1));

  BOOST_CHECK_CLOSE(residual->residual[eBoundLoc0], recoLoc0 - measLoc0, 1e-9);
  BOOST_CHECK_CLOSE(residual->covariance(eBoundLoc0, eBoundLoc0),
                    measVar0 + recoVar0, 1e-9);

  BOOST_CHECK_EQUAL(residual->residual[eBoundLoc1], 0.);
  BOOST_CHECK_EQUAL(residual->covariance(eBoundLoc1, eBoundLoc1), 0.);
}

BOOST_AUTO_TEST_CASE(NoMeasurement) {
  Fixture fixture;
  auto state = fixture.container.trackStateContainer().makeTrackState();
  state.setReferenceSurface(fixture.surface);

  const ConstTrackContainer constContainer(
      std::make_shared<ConstVectorTrackContainer>(fixture.trackBackend),
      std::make_shared<ConstVectorMultiTrajectory>(fixture.trajectory));
  const ConstTrackStateProxy constState =
      constContainer.trackStateContainer().getTrackState(state.index());

  BOOST_CHECK(!measurementResidual(constState, makeParameters(fixture.surface),
                                   TrackParameterType::Predicted)
                   .has_value());
}

BOOST_AUTO_TEST_SUITE_END()
