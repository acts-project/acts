// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

template <typename stepper_t>
class StepperBase
{

public:
  struct State
  {
  };

  /// Jacobian and State defintions
  using Jacobian         = ActsMatrixD<5, 5>;
  using BoundState       = std::tuple<BoundParameters, Jacobian, double>;
  using CurvilinearState = std::tuple<CurvilinearParameters, Jacobian, double>;

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D
  getField(State& state, const Vector3D& pos) const
  {
    return (static_cast<stepper_t*>(this))->getField(state, pos);
  }

  /// Global particle position accessor
  Vector3D
  position(const State& state) const
  {
    return (static_cast<stepper_t*>(this))->position(state);
  }

  /// Momentum direction accessor
  Vector3D
  direction(const State& state) const
  {
    return (static_cast<stepper_t*>(this))->direction(state);
  }

  /// Actual momentum accessor
  double
  momentum(const State& state) const
  {
    return (static_cast<stepper_t*>(this))->momentum(state);
  }

  /// Charge access
  double
  charge(const State& state) const
  {
    return (static_cast<stepper_t*>(this))->charge(state);
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
    return (static_cast<stepper_t*>(this))->surfaceReached(state, surface);
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the surface and creates a bound state. It does not check
  /// if the transported state is at the surface, this needs to
  /// be guaranteed by the propagator
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] reinitialize Boolean flag whether reinitialization is needed,
  /// i.e. if this is an intermediate state of a larger propagation
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  BoundState
  boundState(State&         state,
             const Surface& surface,
             bool           reinitialize = true) const
  {
    return (static_cast<stepper_t*>(this))
        ->boundState(state, surface, reinitialize);
  }

  /// Create and return a curvilinear state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the current position and creates a curvilinear state.
  ///
  /// @param [in] state State that will be presented as @c CurvilinearState
  /// @param [in] reinitialize Boolean flag whether reinitialization is needed
  /// i.e. if this is an intermediate state of a larger propagation
  ///
  /// @return A curvilinear state:
  ///   - the curvilinear parameters at given position
  ///   - the stepweise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  CurvilinearState
  curvilinearState(State& state, bool reinitialize = true) const
  {
    return (static_cast<stepper_t*>(this))
        ->curvilinearState(state, reinitialize);
  }

  /// Method to update a stepper state to the some parameters
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] pars Parameters that will be written into @p state
  void
  update(State& state, const BoundParameters& pars) const
  {
    (static_cast<stepper_t*>(this))->update(state, pars);
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
    (static_cast<stepper_t*>(this))->convert(state, uposition, udirection, up);
  }

  /// Return a corrector
  template <typename corrector_t>
  corrector_t
  corrector(State& state) const
  {
    return (static_cast<stepper_t*>(this))->corrector(state);
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] reinitialize is a flag to steer whether the state should be
  /// reinitialized at the new position
  ///
  /// @return the full transport jacobian
  void
  covarianceTransport(State& state, bool reinitialize = false) const
  {
    (static_cast<stepper_t*>(this))->covarianceTransport(state, reinitialize);
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current position,
  /// or direction of the state
  ///
  /// @tparam surface_t the Surface type
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] surface is the surface to which the covariance is forwarded to
  /// @param [in] reinitialize is a flag to steer whether the state should be
  /// reinitialized at the new position
  /// @note no check is done if the position is actually on the surface
  template <typename surface_t>
  void
  covarianceTransport(State&           state,
                      const surface_t& surface,
                      bool             reinitialize = true) const
  {
    (static_cast<stepper_t*>(this))
        ->covarianceTransport(state, surface, reinitialize);
  }

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param [in,out] state is the propagation state associated with the track
  /// parameters that are being propagated.
  ///
  ///                      the state contains the desired step size.
  ///                      It can be negative during backwards track
  ///                      propagation,
  ///                      and since we're using an adaptive algorithm, it can
  ///                      be modified by the stepper class during propagation.
  template <typename propagator_state_t>
  double
  step(propagator_state_t& state) const
  {
    return (static_cast<stepper_t*>(this))->convert(state);
  }
};
}  // namespace Acts
