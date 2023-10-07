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

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/Surface.hpp"
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
class StraightLineStepper {
 public:
  using Jacobian = BoundMatrix;
  using Covariance = BoundSymMatrix;
  using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<CurvilinearTrackParameters, Jacobian, double>;
  using BField = NullBField;

  /// State for track parameter propagation
  ///
  struct State {
    State() = delete;

    /// Constructor from the initial bound track parameters
    ///
    /// @tparam charge_t Type of the bound parameter charge
    ///
    /// @param [in] gctx is the context object for the geometry
    /// @param [in] mctx is the context object for the magnetic field
    /// @param [in] par The track parameters at start
    /// @param [in] ndir The navigation direction w.r.t momentum
    /// @param [in] ssize is the maximum step size
    /// @param [in] stolerance is the stepping tolerance
    ///
    /// @note the covariance matrix is copied when needed
    template <typename charge_t>
    explicit State(const GeometryContext& gctx,
                   const MagneticFieldContext& mctx,
                   const SingleBoundTrackParameters<charge_t>& par,
                   NavigationDirection ndir = NavigationDirection::Forward,
                   double ssize = std::numeric_limits<double>::max(),
                   double stolerance = s_onSurfaceTolerance)
        : q(par.charge()),
          navDir(ndir),
          stepSize(ndir * std::abs(ssize)),
          tolerance(stolerance),
          geoContext(gctx) {
      (void)mctx;
      pars.template segment<3>(eFreePos0) = par.position(gctx);
      pars.template segment<3>(eFreeDir0) = par.unitDirection();
      pars[eFreeTime] = par.time();
      pars[eFreeQOverP] = par.parameters()[eBoundQOverP];
      if (par.covariance()) {
        // Get the reference surface for navigation
        const auto& surface = par.referenceSurface();
        // set the covariance transport flag to true and copy
        covTransport = true;
        cov = BoundSymMatrix(*par.covariance());
        jacToGlobal = surface.boundToFreeJacobian(gctx, par.parameters());
      }
    }

    /// Jacobian from local to the global frame
    BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

    /// Pure transport jacobian part from runge kutta integration
    FreeMatrix jacTransport = FreeMatrix::Identity();

    /// The full jacobian of the transport entire transport
    Jacobian jacobian = Jacobian::Identity();

    /// The propagation derivative
    FreeVector derivative = FreeVector::Zero();

    /// Internal free vector parameters
    FreeVector pars = FreeVector::Zero();

    /// The charge as the free vector can be 1/p or q/p
    double q = 1.;

    /// Boolean to indiciate if you need covariance transport
    bool covTransport = false;
    Covariance cov = Covariance::Zero();

    /// Navigation direction, this is needed for searching
    NavigationDirection navDir;

    /// accummulated path length state
    double pathAccumulated = 0.;

    /// adaptive step size of the runge-kutta integration
    ConstrainedStep stepSize;

    // Previous step size for overstep estimation (ignored for SL stepper)
    double previousStepSize = 0.;

    /// The tolerance for the stepping
    double tolerance = s_onSurfaceTolerance;

    // Cache the geometry context of this propagation
    std::reference_wrapper<const GeometryContext> geoContext;
  };

  /// Always use the same propagation state type, independently of the initial
  /// track parameter type and of the target surface
  using state_type = State;

  StraightLineStepper() = default;

  template <typename charge_t>
  State makeState(std::reference_wrapper<const GeometryContext> gctx,
                  std::reference_wrapper<const MagneticFieldContext> mctx,
                  const SingleBoundTrackParameters<charge_t>& par,
                  NavigationDirection ndir = NavigationDirection::Forward,
                  double ssize = std::numeric_limits<double>::max(),
                  double stolerance = s_onSurfaceTolerance) const {
    return State{gctx, mctx, par, ndir, ssize, stolerance};
  }

  /// @brief Resets the state
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] boundParams Parameters in bound parametrisation
  /// @param [in] cov Covariance matrix
  /// @param [in] surface The reset @c State will be on this surface
  /// @param [in] navDir Navigation direction
  /// @param [in] stepSize Step size
  void resetState(
      State& state, const BoundVector& boundParams, const BoundSymMatrix& cov,
      const Surface& surface,
      const NavigationDirection navDir = NavigationDirection::Forward,
      const double stepSize = std::numeric_limits<double>::max()) const;

  /// Get the field for the stepping, this gives back a zero field
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Result<Vector3> getField(State& state, const Vector3& pos) const {
    (void)state;
    (void)pos;
    // get the field from the cell
    return Result<Vector3>::success({0., 0., 0.});
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 position(const State& state) const {
    return state.pars.template segment<3>(eFreePos0);
  }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 direction(const State& state) const {
    return state.pars.template segment<3>(eFreeDir0);
  }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double momentum(const State& state) const {
    return std::abs((state.q == 0. ? 1. : state.q) / state.pars[eFreeQOverP]);
  }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double charge(const State& state) const { return state.q; }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double time(const State& state) const { return state.pars[eFreeTime]; }

  /// Overstep limit
  ///
  /// @param state The stepping state (thread-local cache)
  double overstepLimit(const State& state) const {
    (void)state;
    return s_onSurfaceTolerance;
  }

  /// Update surface status
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] surface The surface provided
  /// @param [in] bcheck The boundary check for this status update
  /// @param [in] logger A logger instance
  Intersection3D::Status updateSurfaceStatus(
      State& state, const Surface& surface, const BoundaryCheck& bcheck,
      const Logger& logger = getDummyLogger()) const {
    return detail::updateSingleSurfaceStatus<StraightLineStepper>(
        *this, state, surface, bcheck, logger);
  }

  /// Update step size
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param oIntersection [in] The ObjectIntersection to layer, boundary, etc
  /// @param release [in] boolean to trigger step size release
  template <typename object_intersection_t>
  void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      bool release = true) const {
    detail::updateSingleStepSize<StraightLineStepper>(state, oIntersection,
                                                      release);
  }

  /// Set Step size - explicitely with a double
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The step size value
  /// @param stype [in] The step size type to be set
  /// @param release [in] Do we release the step size?
  void setStepSize(State& state, double stepSize,
                   ConstrainedStep::Type stype = ConstrainedStep::actor,
                   bool release = true) const {
    state.previousStepSize = state.stepSize.value();
    state.stepSize.update(stepSize, stype, release);
  }

  /// Get the step size
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @param stype [in] The step size type to be returned
  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return state.stepSize.value(stype);
  }

  /// Release the Step size
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  void releaseStepSize(State& state) const {
    state.stepSize.release(ConstrainedStep::actor);
  }

  /// Output the Step Size - single component
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  std::string outputStepSize(const State& state) const {
    return state.stepSize.toString();
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief It does not check if the transported state is at the surface, this
  /// needs to be guaranteed by the propagator
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] transportCov Flag steering covariance transport
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  Result<BoundState> boundState(
      State& state, const Surface& surface, bool transportCov = true,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

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
  /// @param [in] freeParams Free parameters that will be written into @p state
  /// @param [in] boundParams Corresponding bound parameters used to update jacToGlobal in @p state
  /// @param [in] covariance Covariance that willl be written into @p state
  /// @param [in] surface The surface used to update the jacToGlobal
  void update(State& state, const FreeVector& freeParams,
              const BoundVector& boundParams, const Covariance& covariance,
              const Surface& surface) const;

  /// Method to update momentum, direction and p
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] up the updated momentum value
  /// @param [in] time the updated time value
  void update(State& state, const Vector3& uposition, const Vector3& udirection,
              double up, double time) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @param [in,out] state State of the stepper
  void transportCovarianceToCurvilinear(State& state) const;

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
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  ///
  void transportCovarianceToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

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
  template <typename propagator_state_t, typename navigator_t>
  Result<double> step(propagator_state_t& state,
                      const navigator_t& /*navigator*/) const {
    // use the adjusted step size
    const auto h = state.stepping.stepSize.value();
    const double p = momentum(state.stepping);
    // time propagates along distance as 1/b = sqrt(1 + m²/p²)
    const auto dtds = std::hypot(1., state.options.mass / p);
    // Update the track parameters according to the equations of motion
    Vector3 dir = direction(state.stepping);
    state.stepping.pars.template segment<3>(eFreePos0) += h * dir;
    state.stepping.pars[eFreeTime] += h * dtds;
    // Propagate the jacobian
    if (state.stepping.covTransport) {
      // The step transport matrix in global coordinates
      FreeMatrix D = FreeMatrix::Identity();
      D.block<3, 3>(0, 4) = ActsSymMatrix<3>::Identity() * h;
      // Extend the calculation by the time propagation
      // Evaluate dt/dlambda
      D(3, 7) = h * state.options.mass * state.options.mass *
                state.stepping.pars[eFreeQOverP] / dtds;
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
