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
    typedef BoundParameters type;
  };

  // ...unless type S is int, in which case it maps to Curvilinear parameters
  template <typename T>
  struct s<T, int>
  {
    typedef CurvilinearParameters type;
  };

public:
  typedef detail::ConstrainedStep cstep;

  /// State for track parameter propagation
  ///
  struct State
  {
    /// Constructor from the initial track parameters
    /// @param[in] par The track parameters at start
    /// @param[in] dir is the navigation direction
    /// @param[in] ssize is the (absolute) maximum step size
    template <typename T>
    explicit State(const T&            par,
                   NavigationDirection ndir = forward,
                   double ssize = std::numeric_limits<double>::max())
      : pos(par.position())
      , dir(par.momentum().normalized())
      , p(par.momentum().norm())
      , q(par.charge())
      , navDir(ndir)
      , pathAccumulated(0.)
      , stepSize(ssize)
    {
    }

    /// Global particle position accessor
    Vector3D
    position() const
    {
      return pos;
    }

    /// Momentum direction accessor
    Vector3D
    direction() const
    {
      return dir;
    }

    /// Momentum accessor
    Vector3D
    momentum() const
    {
      return p * dir;
    }

    /// Charge access
    double
    charge() const
    {
      return q;
    }

    /// Return a corrector
    VoidCorrector
    corrector() const
    {
      return VoidCorrector();
    }

    /// Method to update momentum, direction and p
    ///
    /// @param uposition the updated position
    /// @param udirection the updated direction
    /// @param p the updated momentum value
    void
    update(const Vector3D& uposition, const Vector3D& udirection, double up)
    {
      pos = uposition;
      dir = udirection;
      p   = up;
    }

    /// Global particle position
    Vector3D pos = Vector3D(0, 0, 0);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1, 0, 0);

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
  template <typename T, typename S = int>
  using state_type = State;

  /// Intermediate track parameters are always in curvilinear parametrization
  template <typename T>
  using step_parameter_type = CurvilinearParameters;

  /// Return parameter types depend on the propagation mode:
  /// - when propagating to a surface we return BoundParameters
  /// - otherwise CurvilinearParameters
  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;

  /// Constructor
  StraightLineStepper() = default;

  /// Convert the propagation state (global) to curvilinear parameters
  /// @param state The stepper state
  /// @return curvilinear parameters
  static CurvilinearParameters
  convert(State& state)
  {
    // return the parameters
    return CurvilinearParameters(
        nullptr, state.pos, state.p * state.dir, state.q);
  }

  /// Convert the propagation state to track parameters at a certain surface
  ///
  /// @tparam S The surface type
  ///
  /// @param [in] state Propagation state used
  /// @param [in] surface Destination surface to which the conversion is done
  ///
  /// @return are parameters bound to the target surface
  template <typename S>
  static BoundParameters
  convert(State& state, const S& surface)
  {
    // return the bound parameters
    return BoundParameters(
        nullptr, state.pos, state.p * state.dir, state.q, surface);
  }

  /// Perform a straight line propagation step
  ///
  /// @param[in,out] state is the propagation state associated with the track
  ///                parameters that are being propagated.
  ///                The state contains the desired step size,
  ///                it can be negative during backwards track propagation,
  ///                and since we're using an adaptive algorithm, it can
  ///                be modified by the stepper class during propagation.
  ///
  /// @return the step size taken
  double
  step(State& state) const
  {
    // use the adjusted step size
    const double h = state.stepSize;
    // Update the track parameters according to the equations of motion
    state.pos += h * state.dir;
    // state the path length
    state.pathAccumulated += h;
    // return h
    return h;
  }

  /// Get the field for the stepping, this gives back a zero field
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D
  getField(State&, const Vector3D&) const
  {
    // get the field from the cell
    return Vector3D(0., 0., 0.);
  }
};

}  // namespace Acts
