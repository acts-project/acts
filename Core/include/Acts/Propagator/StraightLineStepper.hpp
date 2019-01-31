// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// @brief straight line stepper based on Surface intersection
///
/// The straight line stepper is a simple navigation stepper
/// to be used to navigate through the tracking geometry. It can be
/// used for simple material mapping, navigation validation
class StraightLineStepper
{

private:
  // This struct is a meta-function which normally maps to BoundParameters...
  template <typename T, typename S>
  struct s
  {
    using type = BoundParameters;
  };

  // ...unless type S is int, in which case it maps to Curvilinear parameters
  template <typename T>
  struct s<T, int>
  {
    using type = CurvilinearParameters;
  };

public:
  using cstep = detail::ConstrainedStep;

  /// State for track parameter propagation
  ///
  struct State
  {

    /// Constructor from the initial track parameters
    ///
    /// @tparam parameters_t the Type of the track parameters
    ///
    /// @param [in] par The track parameters at start
    /// @param [in] dir is the navigation direction
    /// @param [in] ssize is the (absolute) maximum step size
    template <typename parameters_t>
    explicit State(const parameters_t& par,
                   NavigationDirection ndir = forward,
                   double ssize = std::numeric_limits<double>::max())
      : pos(par.position())
      , dir(par.momentum().normalized())
      , p(par.momentum().norm())
      , q(par.charge())
      , navDir(ndir)
      , stepSize(ssize)
    {
    }

    /// Boolean to indiciate if you need covariance transport
    bool              covTransport = false;
    ActsSymMatrixD<5> cov          = ActsSymMatrixD<5>::Zero();

    /// Global particle position
    Vector3D pos = Vector3D(0., 0., 0.);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1., 0., 0.);

    /// Momentum
    double p = 0.;

    /// Save the charge: neutral as default for SL stepper
    double q = 0.;

    /// Navigation direction, this is needed for searching
    NavigationDirection navDir;

    /// accummulated path length state
    double pathAccumulated = 0.;

    /// adaptive step size of the runge-kutta integration
    cstep stepSize = std::numeric_limits<double>::max();
  };

  /// Global particle position accessor
  Vector3D
  position(const State& state) const
  {
    return state.pos;
  }

  /// Momentum direction accessor
  Vector3D
  direction(const State& state) const
  {
    return state.dir;
  }

  /// Momentum accessor
  double
  momentum(const State& state) const
  {
    return state.p;
  }

  /// Charge access
  double
  charge(const State& state) const
  {
    return state.q;
  }

  /// Return a corrector
  static VoidIntersectionCorrector
  corrector(State& /*state*/)
  {
    return VoidIntersectionCorrector();
  }

  /// Method to update momentum, direction and p
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] up the updated momentum value
  static void
  update(State&          state,
         const Vector3D& uposition,
         const Vector3D& udirection,
         double          up)
  {
    state.pos = uposition;
    state.dir = udirection;
    state.p   = up;
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] reinitialize is a flag to steer whether the
  ///        state should be reinitialized at the new
  ///        position
  ///
  /// @return the full transport jacobian
  static const ActsMatrixD<5, 5>
  covarianceTransport(State& /*state*/, bool /*reinitialize = false*/)
  {
    return ActsMatrixD<5, 5>::Identity();
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @tparam surface_t the surface type - ignored here
  ///
  /// @param [in,out] state The stepper state
  /// @param [in] surface is the surface to which the covariance is
  ///        forwarded to
  /// @param [in] reinitialize is a flag to steer whether the
  ///        state should be reinitialized at the new
  ///        position
  /// @note no check is done if the position is actually on the surface
  ///
  /// @return the full transport jacobian
  template <typename surface_t>
  static const ActsMatrixD<5, 5>
  covarianceTransport(State& /*unused*/,
                      const surface_t& /*surface*/,
                      bool /*reinitialize = false*/)
  {
    return ActsMatrixD<5, 5>::Identity();
  }

  /// Always use the same propagation state type, independently of the initial
  /// track parameter type and of the target surface
  template <typename parameters_t, typename surface_t = int>
  using state_type = State;

  /// Intermediate track parameters are always in curvilinear parametrization
  template <typename parameters_t>
  using step_parameter_type = CurvilinearParameters;

  /// Return parameter types depend on the propagation mode:
  /// - when propagating to a surface we return BoundParameters
  /// - otherwise CurvilinearParameters
  template <typename parameters_t, typename surface_t = int>
  using return_parameter_type = typename s<parameters_t, surface_t>::type;

  /// Constructor
  StraightLineStepper() = default;

  /// Convert the propagation state (global) to curvilinear parameters
  ///
  /// @tparam result_t Type of the propagator result to be filled
  ///
  /// @param [in] state The stepper state
  /// @param [in,out] result The result object from the propagator
  template <typename result_t>
  void
  convert(State& state, result_t& result) const
  {
    // Fill the end parameters
    result.endParameters = std::make_unique<const CurvilinearParameters>(
        nullptr, state.pos, state.p * state.dir, state.q);
  }

  /// Convert the propagation state to track parameters at a certain surface
  ///
  /// @tparam result_t Type of the propagator result to be filled
  /// @tparam surface_t Type of the surface
  ///
  /// @param [in,out] state Propagation state used
  /// @param [in,out] result Result object from the propagator
  /// @param [in] surface Destination surface to which the conversion is done
  template <typename result_t, typename surface_t>
  void
  convert(State& state, result_t& result, const surface_t& surface) const
  {
    // Fill the end parameters
    result.endParameters
        = std::make_unique<const BoundParameters>(nullptr,
                                                  state.pos,
                                                  state.p * state.dir,
                                                  state.q,
                                                  surface.getSharedPtr());
  }

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
  double
  step(propagator_state_t& state) const
  {
    // use the adjusted step size
    const double h = state.stepping.stepSize;
    // Update the track parameters according to the equations of motion
    state.stepping.pos += h * state.stepping.dir;
    // state the path length
    state.stepping.pathAccumulated += h;
    // return h
    return h;
  }

  /// Get the field for the stepping, this gives back a zero field
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D
  getField(State& /*state*/, const Vector3D& /*pos*/) const
  {
    // get the field from the cell
    return Vector3D(0., 0., 0.);
  }
};

}  // namespace Acts
