// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/ConstrainedStepControl.hpp"
#include "Acts/Propagator/StepperState.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <functional>

namespace Acts {

/// @brief straight line stepper based on Surface intersection
///
/// The straight line stepper is a simple navigation stepper
/// to be used to navigate through the tracking geometry. It can be
/// used for simple material mapping, navigation validation
///
/// A nested State object guarantees for the thread-local caching of the
/// transport parameters. This state object is provided by the Propagator
/// and the stepper is used to interpret the non-trivial internal parameters,
/// e.g. the current global position is gathered through:
///
///   auto position = stepper.position(state);
///
class StraightLineStepper {
 public:
  using Jacobian = BoundMatrix;
  using Covariance = BoundSymMatrix;
  using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<CurvilinearTrackParameters, Jacobian, double>;
  using BField = NullBField;

  /// Always use the same propagation state type, independently of the initial
  /// track parameter type and of the target surface
  using State = StepperState<StraightLineStepper>;

  /// The control mechanism for the constrained step size
  using StepControl = ConstrainedStepControl<StraightLineStepper>;
  StepControl stepControl;

  StraightLineStepper() = default;

  /// @brief Resets the state
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] boundParams Parameters in bound parametrisation
  /// @param [in] freeParams Parameters in free parametrisation
  /// @param [in] cov Covariance matrix
  /// @param [in] navDir Navigation direction
  /// @param [in] stepSize Step size
  void resetState(
      State& state, const BoundVector& boundParams, const BoundSymMatrix& cov,
      const Surface& surface, const NavigationDirection navDir = forward,
      const double stepSize = std::numeric_limits<double>::max()) const;

  /// Get the field for the stepping, this gives back a zero field
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D getField(State& /*state*/, const Vector3D& /*pos*/) const {
    // get the field from the cell
    return Vector3D(0., 0., 0.);
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3D position(const State& state) const {
    return state.pars.template segment<3>(eFreePos0);
  }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3D direction(const State& state) const {
    return state.pars.template segment<3>(eFreeDir0);
  }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double momentum(const State& state) const {
    return std::abs(1. / state.pars[eFreeQOverP]);
  }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double charge(const State& state) const { return state.q; }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double time(const State& state) const { return state.pars[eFreeTime]; }

  /// Access to the current geometry context
  ///
  /// @param state [in] The stepping state (thread-local cache)
  const GeometryContext& geometryContext(const State& state) const {
    return state.geoContext;
  }

  /// Access to the navigation direction
  ///
  /// @param state [in] The stepping state (thread-local cache)
  NavigationDirection steppingDirection(const State& state) const {
    return state.navDir;
  }

  /// Access to the navigation direction
  ///
  /// @param state [in, out] The stepping state (thread-local cache)
  /// @param sdir [in] stepping direction
  void setSteppingDirection(State& state, NavigationDirection sdir) const {
    state.navDir = sdir;
  }

  /// Access to the stepping tolerance
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double steppingTolerance(const State& state) const { return state.tolerance; }

  /// Access to the accumulated path
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double accumulatedPath(const State& state) const {
    return state.pathAccumulated;
  }

  /// Reset to the accumulated path
  ///
  /// @param state [in, out] The stepping state (thread-local cache)
  void resetAccumulatedPath(State& state) const {
    state.pathAccumulated = 0;
    ;
  }

  /// Indicate if the covariance has to be transported
  ///
  /// @param state [in] The stepping state (thread-local cache)
  bool transportCovariance(const State& state) const {
    return state.covTransport;
  }

  /// Overstep limit
  ///
  /// @param state The stepping state (thread-local cache)
  double overstepLimit(const State& /*state*/) const {
    return s_onSurfaceTolerance;
  }

  /// Update surface status
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param surface [in] The surface provided
  /// @param bcheck [in] The boundary check for this status update
  Intersection3D::Status updateSurfaceStatus(
      State& state, const Surface& surface, const BoundaryCheck& bcheck) const {
    return detail::updateSingleSurfaceStatus<StraightLineStepper>(
        *this, state, surface, bcheck);
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief It does not check if the transported state is at the surface, this
  /// needs to be guaranteed by the propagator
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  BoundState boundState(State& state, const Surface& surface,
                        bool transportCov = true) const;

  /// Create and return a curvilinear state at the current position
  ///
  /// @brief This creates a curvilinear state.
  ///
  /// @param [in] state State that will be presented as @c CurvilinearState
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A curvilinear state:
  ///   - the curvilinear parameters at given position
  ///   - the stepweise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  CurvilinearState curvilinearState(State& state,
                                    bool transportCov = true) const;

  /// Method to update a stepper state to the some parameters
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] pars Parameters that will be written into @p state
  void update(State& state, const FreeVector& parameters,
              const Covariance& covariance) const;

  /// Method to update momentum, direction and p
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] up the updated momentum value
  /// @param [in] time the updated time value
  void update(State& state, const Vector3D& uposition,
              const Vector3D& udirection, double up, double time) const;

  /// @brief Convenience method for better readability
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] change The delta that may be applied to it
  ///
  void updateBoundVariance(State& state, BoundIndices bIndex,
                           double delta) const;

  /// @brief Convenience method for better readability
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] cov the bound covariance matrix to be updated
  ///
  void updateBoundCovariance(State& state, Covariance cov) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @param [in,out] state State of the stepper
  void covarianceTransport(State& state) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @tparam surface_t the surface type - ignored here
  ///
  /// @param [in,out] state The stepper state
  /// @param [in] surface is the surface to which the covariance is
  ///        forwarded to
  /// @note no check is done if the position is actually on the surface
  ///
  void covarianceTransport(State& state, const Surface& surface) const;

  /// Perform a straight line propagation step
  ///
  /// @param [in,out] state is the propagation state associated with the track
  /// parameters that are being propagated.
  ///                The state contains the desired step size,
  ///                it can be negative during backwards track propagation,
  ///                and since we're using an adaptive algorithm, it can
  ///                be modified by the stepper class during propagation.
  ///
  /// @return the step size taken
  template <typename propagator_state_t>
  Result<double> step(propagator_state_t& state) const {
    // use the adjusted step size
    const auto h = state.stepping.stepSize;
    const double p = momentum(state.stepping);
    // time propagates along distance as 1/b = sqrt(1 + m²/p²)
    const auto dtds = std::hypot(1., state.options.mass / p);
    // Update the track parameters according to the equations of motion
    Vector3D dir = direction(state.stepping);
    state.stepping.pars.template segment<3>(eFreePos0) += h * dir;
    state.stepping.pars[eFreeTime] += h * dtds;
    // Propagate the jacobian
    if (state.stepping.covTransport) {
      // The step transport matrix in global coordinates
      FreeMatrix D = FreeMatrix::Identity();
      D.block<3, 3>(0, 4) = ActsSymMatrixD<3>::Identity() * h;
      // Extend the calculation by the time propagation
      // Evaluate dt/dlambda
      D(3, 7) = h * state.options.mass * state.options.mass *
                (state.stepping.q == 0. ? 1. : state.stepping.q) / (p * dtds);
      // Set the derivative factor the time
      state.stepping.derivative(3) = dtds;
      // Update jacobian and derivative
      state.stepping.jacTransport = D * state.stepping.jacTransport;
      state.stepping.derivative.template head<3>() = dir;
    }
    // state the path length
    state.stepping.pathAccumulated += h;

    // return h
    return h;
  }
};

}  // namespace Acts
