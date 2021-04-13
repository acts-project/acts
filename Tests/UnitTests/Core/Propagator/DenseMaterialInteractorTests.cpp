// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Propagator/DenseMaterialInteractor.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

static constexpr auto eps = 3 * std::numeric_limits<double>::epsilon();

/// Simplified stepper state
struct StepperState {
  BoundSymMatrix cov = BoundSymMatrix::Identity();
  ActsScalar pathAccumulated;
  FreeMatrix jacTransport;

  ActsScalar charge;
  ActsScalar momentum;
  ActsScalar time;
  Vector3 direction;
  Vector3 position;
};

/// Simplified propagator state
struct State {
  struct {
    TrackingVolume const* currentVolume;
    bool targetReached;
  } navigation;

  StepperState stepping;

  struct {
    int absPdgCode;
    ActsScalar mass;
  } options;
};

/// Simplified stepper
struct Stepper {
  void covarianceTransport(StepperState& state) const {
    state.jacTransport = FreeMatrix::Identity();
  }
  ActsScalar charge(StepperState& state) const { return state.charge; }
  ActsScalar momentum(StepperState& state) const { return state.momentum; }
  ActsScalar time(StepperState& state) const { return state.time; }
  Vector3 direction(StepperState& state) const { return state.direction; }
  Vector3 position(StepperState& state) const { return state.position; }
};

BOOST_AUTO_TEST_CASE(test_material_collector) {
  // Build the required objects
  DenseMaterialInteractor dmi;
  DenseMaterialInteractor::Result result;
  State state;
  Stepper stepper;

  /// Case 1: No if-condition is fulfilled
  /// a) Due to missing volume
  // Construct the state
  state.navigation.currentVolume = nullptr;
  result.currentVolume = nullptr;
  state.stepping.pathAccumulated = 1.;
  state.stepping.position = Vector3(2., 3., 4.);
  state.stepping.jacTransport = FreeSymMatrix::Identity() * 2.;
  state.navigation.targetReached = false;
  // Call the actor
  dmi(state, stepper, result);
  // Check the resuts
  BOOST_CHECK_EQUAL(result.s, 0.);
  BOOST_CHECK_EQUAL(result.qop0, 0.);
  BOOST_CHECK_EQUAL(result.t0, 0.);
  BOOST_CHECK_EQUAL(result.currentVolume, state.navigation.currentVolume);
  BOOST_CHECK_EQUAL(result.sAccumulated, state.stepping.pathAccumulated);
  BOOST_CHECK_EQUAL(result.sInX0, 0.);
  BOOST_CHECK_EQUAL(result.dir0, Vector3::Zero());
  BOOST_CHECK_EQUAL(result.previousPos, state.stepping.position);
  BOOST_CHECK_EQUAL(state.stepping.cov, BoundSymMatrix::Identity());

  /// b) Due to missing material
  // Build the volumes without material
  Transform3 transform = Transform3::Identity();
  std::shared_ptr<const VolumeBounds> volumeBounds =
      std::make_shared<const CuboidVolumeBounds>(1_m, 1_m, 1_m);
  TrackingVolume const* trackingVolume =
      TrackingVolume::create(transform, volumeBounds).get();
  state.navigation.currentVolume = trackingVolume;
  result.currentVolume = trackingVolume;
  // Change the state a bit
  state.stepping.pathAccumulated *= 2.;
  state.stepping.position *= 2.;
  // Call the actor
  dmi(state, stepper, result);
  // Check the resuts
  BOOST_CHECK_EQUAL(result.s, 0.);
  BOOST_CHECK_EQUAL(result.qop0, 0.);
  BOOST_CHECK_EQUAL(result.t0, 0.);
  BOOST_CHECK_EQUAL(result.currentVolume, state.navigation.currentVolume);
  BOOST_CHECK_EQUAL(result.sAccumulated, state.stepping.pathAccumulated);
  BOOST_CHECK_EQUAL(result.sInX0, 0.);
  BOOST_CHECK_EQUAL(result.dir0, Vector3::Zero());
  BOOST_CHECK_EQUAL(result.previousPos, state.stepping.position);
  BOOST_CHECK_EQUAL(state.stepping.cov, BoundSymMatrix::Identity());

  /// Case 2: The material is recorded
  /// a) Plain recording without adding to the covariance matrix
  // Build a volume with material
  Material material = makeSilicon();
  std::shared_ptr<const IVolumeMaterial> volumeMaterial =
      std::make_shared<const HomogeneousVolumeMaterial>(material);
  TrackingVolume const* materialTrackingVolume =
      TrackingVolume::create(transform, volumeBounds, volumeMaterial).get();
  state.navigation.currentVolume = materialTrackingVolume;
  result.currentVolume = materialTrackingVolume;
  // Change the state a bit
  state.stepping.pathAccumulated *= 2.;
  state.stepping.position *= 2.;
  // Call the actor
  dmi(state, stepper, result);
  // Check the resuts
  CHECK_CLOSE_ABS(result.s, 0.5 * state.stepping.pathAccumulated, eps);
  BOOST_CHECK_EQUAL(result.currentVolume, state.navigation.currentVolume);
  BOOST_CHECK_EQUAL(result.sAccumulated, state.stepping.pathAccumulated);
  CHECK_CLOSE_ABS(result.sInX0,
                  0.5 * state.stepping.pathAccumulated / material.X0(), eps);
  BOOST_CHECK_EQUAL(result.previousPos, state.stepping.position);
  BOOST_CHECK_EQUAL(state.stepping.cov, BoundSymMatrix::Identity());

  /// b) We have some material and repeat it - no covariance modification
  // Go a bit further
  state.stepping.pathAccumulated *= 2.;
  state.stepping.position *= 2.;
  // Call the actor
  dmi(state, stepper, result);
  // Check the resuts
  CHECK_CLOSE_ABS(result.s, 0.75 * state.stepping.pathAccumulated, eps);
  BOOST_CHECK_EQUAL(result.currentVolume, state.navigation.currentVolume);
  BOOST_CHECK_EQUAL(result.sAccumulated, state.stepping.pathAccumulated);
  CHECK_CLOSE_ABS(result.sInX0,
                  0.75 * state.stepping.pathAccumulated / material.X0(), eps);
  BOOST_CHECK_EQUAL(result.previousPos, state.stepping.position);
  BOOST_CHECK_EQUAL(state.stepping.cov, BoundSymMatrix::Identity());

  /// Case 3: Initialisation from the state by changing the volume
  // Configure the states
  result.s = 0.;
  result.currentVolume = nullptr;
  state.stepping.charge = 1_e;
  state.stepping.momentum = 42_GeV;
  state.stepping.time = 1337_s;
  state.stepping.direction = Vector3(1., 2.5, 33.33).normalized();
  // Call the actor
  dmi(state, stepper, result);
  // Check the resuts
  BOOST_CHECK_EQUAL(result.s, 0.);
  CHECK_CLOSE_ABS(result.qop0, state.stepping.charge / state.stepping.momentum,
                  eps);
  BOOST_CHECK_EQUAL(result.t0, state.stepping.time);
  BOOST_CHECK_EQUAL(result.sInX0, 0.);
  BOOST_CHECK_EQUAL(result.dir0, state.stepping.direction);
  BOOST_CHECK_EQUAL(result.currentVolume, state.navigation.currentVolume);
  BOOST_CHECK_EQUAL(result.sAccumulated, state.stepping.pathAccumulated);
  BOOST_CHECK_EQUAL(state.stepping.cov, BoundSymMatrix::Identity());

  /// Case 4: Adding the scattering contribution to the covariance
  /// a) Due to changing the volume
  // Configure the states
  result.s = 9001.;
  result.sInX0 = 3.;
  result.currentVolume = nullptr;
  // Call the actor
  dmi(state, stepper, result);
  // Check the resuts
  BOOST_CHECK_EQUAL(result.s, 0.);
  CHECK_CLOSE_ABS(result.qop0, state.stepping.charge / state.stepping.momentum,
                  eps);
  BOOST_CHECK_EQUAL(result.t0, state.stepping.time);
  BOOST_CHECK_EQUAL(result.sInX0, 0.);
  BOOST_CHECK_EQUAL(result.dir0, state.stepping.direction);
  BOOST_CHECK_EQUAL(result.currentVolume, state.navigation.currentVolume);
  BOOST_CHECK_EQUAL(result.sAccumulated, state.stepping.pathAccumulated);
  BOOST_CHECK_NE(state.stepping.cov, BoundSymMatrix::Identity());

  /// b) Due to a triggered covariance transport from a different place
  // Configure the states
  result.s = 9001.;
  result.sInX0 = 3.;
  state.stepping.cov = BoundSymMatrix::Identity();
  state.stepping.jacTransport = FreeMatrix::Identity();
  // Call the actor
  dmi(state, stepper, result);
  // Check the resuts
  BOOST_CHECK_EQUAL(result.s, 0.);
  CHECK_CLOSE_ABS(result.qop0, state.stepping.charge / state.stepping.momentum,
                  eps);
  BOOST_CHECK_EQUAL(result.t0, state.stepping.time);
  BOOST_CHECK_EQUAL(result.sInX0, 0.);
  BOOST_CHECK_EQUAL(result.dir0, state.stepping.direction);
  BOOST_CHECK_EQUAL(result.currentVolume, state.navigation.currentVolume);
  BOOST_CHECK_EQUAL(result.sAccumulated, state.stepping.pathAccumulated);
  BOOST_CHECK_NE(state.stepping.cov, BoundSymMatrix::Identity());

  /// c) Due to target reached
  // Configure the states
  result.s = 9001.;
  result.sInX0 = 3.;
  state.stepping.cov = BoundSymMatrix::Identity();
  state.stepping.jacTransport = FreeMatrix::Identity() * 2.;
  state.navigation.targetReached = true;
  // Call the actor
  dmi(state, stepper, result);
  // Check the resuts
  BOOST_CHECK_EQUAL(result.s, 0.);
  CHECK_CLOSE_ABS(result.qop0, state.stepping.charge / state.stepping.momentum,
                  eps);
  BOOST_CHECK_EQUAL(result.t0, state.stepping.time);
  BOOST_CHECK_EQUAL(result.sInX0, 0.);
  BOOST_CHECK_EQUAL(result.dir0, state.stepping.direction);
  BOOST_CHECK_EQUAL(result.currentVolume, state.navigation.currentVolume);
  BOOST_CHECK_EQUAL(result.sAccumulated, state.stepping.pathAccumulated);
  BOOST_CHECK_NE(state.stepping.cov, BoundSymMatrix::Identity());
}

}  // namespace Test
}  // namespace Acts
