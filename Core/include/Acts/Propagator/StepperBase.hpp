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

  using Jacobian         = ActsMatrixD<5, 5>;
  using BoundState       = std::tuple<BoundParameters, Jacobian, double>;
  using CurvilinearState = std::tuple<CurvilinearParameters, Jacobian, double>;

  template <typename result_t>
  void
  convert(State& state, result_t& result) const
  {
    (static_cast<stepper_t*>(this))->convert(state, result);
  }

  template <typename result_t, typename surface_t>
  void
  convert(State& state, result_t& result, const surface_t& surface) const
  {
    (static_cast<stepper_t*>(this))->convert(state, result);
  }

  Vector3D
  getField(State& state, const Vector3D& pos) const
  {
    return (static_cast<stepper_t*>(this))->getField(state, pos);
  }

  Vector3D
  position(const State& state) const
  {
    return (static_cast<stepper_t*>(this))->position(state);
  }

  Vector3D
  direction(const State& state) const
  {
    return (static_cast<stepper_t*>(this))->direction(state);
  }

  double
  momentum(const State& state) const
  {
    return (static_cast<stepper_t*>(this))->momentum(state);
  }

  double
  charge(const State& state) const
  {
    return (static_cast<stepper_t*>(this))->charge(state);
  }

  bool
  surfaceReached(const State& state, const Surface* surface) const
  {
    return (static_cast<stepper_t*>(this))->surfaceReached(state, surface);
  }

  BoundState
  boundState(State&         state,
             const Surface& surface,
             bool           reinitialize = true) const
  {
    return (static_cast<stepper_t*>(this))
        ->boundState(state, surface, reinitialize);
  }

  CurvilinearState
  curvilinearState(State& state, bool reinitialize = true) const
  {
    return (static_cast<stepper_t*>(this))
        ->curvilinearState(state, reinitialize);
  }

  void
  update(State& state, const BoundParameters& pars) const
  {
    (static_cast<stepper_t*>(this))->update(state, pars);
  }

  void
  update(State&          state,
         const Vector3D& uposition,
         const Vector3D& udirection,
         double          up) const
  {
    (static_cast<stepper_t*>(this))->convert(state, uposition, udirection, up);
  }

  template <typename corrector_t>
  corrector_t
  corrector(State& state) const
  {
    return (static_cast<stepper_t*>(this))->corrector(state);
  }

  void
  covarianceTransport(State& state, bool reinitialize = false) const
  {
    (static_cast<stepper_t*>(this))->covarianceTransport(state, reinitialize);
  }

  template <typename surface_t>
  void
  covarianceTransport(State&           state,
                      const surface_t& surface,
                      bool             reinitialize = true) const
  {
    (static_cast<stepper_t*>(this))
        ->covarianceTransport(state, surface, reinitialize);
  }

  template <typename propagator_state_t>
  double
  step(propagator_state_t& state) const
  {
    return (static_cast<stepper_t*>(this))->convert(state);
  }
};
}  // namespace Acts
