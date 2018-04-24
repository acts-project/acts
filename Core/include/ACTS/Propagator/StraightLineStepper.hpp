// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_STRAIGHTLINE_STEPPER_HPP
#define ACTS_STRAIGHTLINE_STEPPER_HPP

#include <cmath>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/MagneticField/concept/AnyFieldLookup.hpp"
#include "ACTS/Propagator/detail/ConstrainedStep.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

#ifndef SLTEPPER_DEBUG_OUTPUTS
#define SLTEPPER_DEBUG_OUTPUTS
#define SLLOG(cache, message)                                                  \
  if (cache.debug) {                                                           \
    std::stringstream dstream;                                                 \
    dstream << "|->" << std::setw(cache.debugPfxWidth);                        \
    dstream << "StraightLineStepper"                                           \
            << " | ";                                                          \
    dstream << std::setw(cache.debugMsgWidth) << message << '\n';              \
    cache.debugString += dstream.str();                                        \
  }
#endif

namespace Acts {

/// StraightLineStepper
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

  /// Cache for track parameter propagation
  ///
  struct Cache
  {
    /// Constructor from the initial track parameters
    /// @param [in] par The track parameters at start
    ///
    /// @note the covariance matrix is copied when needed
    template <typename T>
    explicit Cache(const T&            par,
                   NavigationDirection ndir = forward,
                   double ssize = std::numeric_limits<double>::max())
      : pos(par.position())
      , dir(par.momentum().normalized())
      , qop(par.charge() / par.momentum().norm())
      , navDir(ndir)
      , accumulatedPath(0.)
      , stepSize(ndir * ssize)
    {
      // Get the reference surface for navigation
      const auto& surface = par.referenceSurface();
      // cache the surface for navigation
      startSurface = &surface;
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

    /// Actual momentum accessor
    Vector3D
    momentum() const
    {
      return (1. / qop) * dir;
    }

    /// Global particle position
    Vector3D pos = Vector3D(0, 0, 0);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1, 0, 0);

    /// Charge-momentum ratio, in natural units
    double qop = 1;

    /// Navigation direction, this is needed for searching
    NavigationDirection navDir;

    /// accummulated path length cache
    double accumulatedPath = 0.;

    /// adaptive step size of the runge-kutta integration
    cstep stepSize = std::numeric_limits<double>::max();

    /// Navigation cache: the start surface
    const Surface* startSurface = nullptr;

    /// Navigation cache: the current surface
    const Surface* currentSurface = nullptr;

    /// Navigation cache: the target surface
    const Surface* targetSurface = nullptr;
    bool           targetReached = false;

    /// Debug output
    /// the string where things are stored (optionally)
    bool        debug       = false;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    size_t debugPfxWidth = 30;
    size_t debugMsgWidth = 50;
  };

  /// Always use the same propagation cache type, independently of the initial
  /// track parameter type and of the target surface
  template <typename T, typename S = int>
  using cache_type = Cache;

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

  /// Convert the propagation cache (global) to curvilinear parameters
  /// @param cache The stepper cache
  /// @return curvilinear parameters
  static CurvilinearParameters
  convert(Cache& cache)
  {
    double charge = cache.qop > 0. ? 1. : -1.;
    // return the parameters
    return CurvilinearParameters(
        nullptr, cache.pos, cache.dir / std::abs(cache.qop), charge);
  }

  /// Convert the propagation cache to track parameters at a certain surface
  ///
  /// @tparam S The surface type
  ///
  /// @param [in] cache Propagation cache used
  /// @param [in] surface Destination surface to which the conversion is done
  ///
  /// @return are parameters bound to the target surface
  template <typename S>
  static BoundParameters
  convert(Cache& cache, const S& surface)
  {
    double charge = cache.qop > 0. ? 1. : -1.;
    // return the bound parameters
    return BoundParameters(
        nullptr, cache.pos, cache.dir / std::abs(cache.qop), charge, surface);
  }

  /// Perform a straight line propagation step
  ///
  /// @param[in,out] cache is the propagation cache associated with the track
  ///                parameters that are being propagated.
  ///                The cache contains the desired step size,
  ///                it can be negative during backwards track propagation,
  ///                and since we're using an adaptive algorithm, it can
  ///                be modified by the stepper class during propagation.
  ///
  /// @return the step size taken
  double
  step(Cache& cache) const
  {
    // use the adjusted step size
    const double h = cache.stepSize;
    // debug output
    SLLOG(cache, "Performing StraightLine step with size " << h);
    // Update the track parameters according to the equations of motion
    cache.pos += h * cache.dir;
    // cache the path length
    cache.accumulatedPath += h;
    // return h
    return h;
  }
};

}  // namespace Acts

#endif  // ACTS_STRAIGHTLINE_STEPPER_HPP
