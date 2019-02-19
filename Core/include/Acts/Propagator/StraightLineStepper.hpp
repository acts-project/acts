// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Context.hpp"
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

  using Jacobian         = ActsMatrixD<5, 5>;
  using BoundState       = std::tuple<BoundParameters, Jacobian, double>;
  using CurvilinearState = std::tuple<CurvilinearParameters, Jacobian, double>;

  /// State for track parameter propagation
  ///
  struct State
  {

    /// Constructor from the initial track parameters
    ///
    /// @tparam parameters_t the Type of the track parameters
    ///
    /// @param [in] ctx is the context object, ingored
    /// @param [in] par The track parameters at start
    /// @param [in] dir is the navigation direction
    /// @param [in] ssize is the (absolute) maximum step size
    template <typename parameters_t>
    explicit State(Context /*ctx*/,
                   const parameters_t& par,
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

  /// Always use the same propagation state type, independently of the initial
  /// track parameter type and of the target surface
  using state_type = State;

  /// Return parameter types depend on the propagation mode:
  /// - when propagating to a surface we return BoundParameters
  /// - otherwise CurvilinearParameters
  template <typename parameters_t, typename surface_t = int>
  using return_parameter_type = typename s<parameters_t, surface_t>::type;

  /// Constructor
  StraightLineStepper() = default;

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

  /// Tests if the state reached a surface
  ///
  /// @param [in] state State that is tests
  /// @param [in] surface Surface that is tested
  ///
  /// @return Boolean statement if surface is reached by state
  bool
  surfaceReached(const State& state, const Surface* surface) const
  {
    return surface->isOnSurface(position(state), direction(state), true);
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief It does not check if the transported state is at the surface, this
  /// needs to be guaranteed by the propagator
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  BoundState
  boundState(State& state, const Surface& surface, bool /*unused*/) const
  {
    // Create the bound parameters
    BoundParameters parameters(nullptr,
                               state.pos,
                               state.p * state.dir,
                               state.q,
                               surface.getSharedPtr());
    // Create the bound state
    BoundState bState{std::move(parameters),
                      ActsMatrixD<5, 5>::Identity(),
                      state.pathAccumulated};
    /// Return the State
    return bState;
  }

  /// Create and return a curvilinear state at the current position
  ///
  /// @brief This creates a curvilinear state.
  ///
  /// @param [in] state State that will be presented as @c CurvilinearState
  ///
  /// @return A curvilinear state:
  ///   - the curvilinear parameters at given position
  ///   - the stepweise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  CurvilinearState
  curvilinearState(State& state, bool /*unused*/) const
  {
    // Create the curvilinear parameters
    CurvilinearParameters parameters(
        nullptr, state.pos, state.p * state.dir, state.q);
    // Create the bound state
    CurvilinearState curvState{std::move(parameters),
                               ActsMatrixD<5, 5>::Identity(),
                               state.pathAccumulated};
    /// Return the State
    return curvState;
  }

  /// Method to update a stepper state to the some parameters
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] pars Parameters that will be written into @p state
  void
  update(State& state, const BoundParameters& pars) const
  {
    const auto& mom = pars.momentum();
    state.pos       = pars.position();
    state.dir       = mom.normalized();
    state.p         = mom.norm();
  }

  /// Method to update momentum, direction and p
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] up the updated momentum value
  void
  update(State&          state,
         const Vector3D& uposition,
         const Vector3D& udirection,
         double          up) const
  {
    state.pos = uposition;
    state.dir = udirection;
    state.p   = up;
  }

  /// Return a corrector
  VoidIntersectionCorrector
  corrector(State& /*state*/) const
  {
    return VoidIntersectionCorrector();
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] reinitialize is a flag to steer whether the
  ///        state should be reinitialized at the new
  ///        position
  void
  covarianceTransport(State& /*state*/, bool /*reinitialize = false*/) const
  {
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @tparam surface_t the surface type - ignored here
  ///
  /// @param [in] ctx The context for thsi call - ingnored here
  /// @param [in,out] state The stepper state
  /// @param [in] surface is the surface to which the covariance is
  ///        forwarded to
  /// @param [in] reinitialize is a flag to steer whether the
  ///        state should be reinitialized at the new
  ///        position
  /// @note no check is done if the position is actually on the surface
  ///
  /// @return the full transport jacobian
  static const ActsMatrixD<5, 5>
  covarianceTransport(Context /*ctx*/,
                      State& /*unused*/,
                      const Surface& /*surface*/,
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
};

}  // namespace Acts
