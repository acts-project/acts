// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <sstream>
#include <string>
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

namespace detail {

  using Cstep = ConstrainedStep;

  /// @brief a struct that allows you to evaluate the intersection
  ///
  /// It actually modifies the position/direction of the intersection
  /// attempt to allow for a potential better estimate when being
  /// in a magnetic field setup
  ///
  /// This is acutally most relevant for a constant magnetic field where
  /// the danger of overstepping is problematic
  struct RelativePathCorrector
  {

    /// Start position where this corrector is created
    Vector3D startPos = Vector3D(0., 0., 0.);
    /// Start direction with which this corrector is created
    Vector3D startDir = Vector3D(0., 0., 0.);
    /// Path length where this corrector is created
    double pathLength = 0.;
    /// Below here do only straight line estimate
    double straightLineStep = 100 * units::_um;
    /// Step modification factor
    double stepModification = 0.5;

    /// Constructor from arguments
    RelativePathCorrector(const Vector3D& spos  = Vector3D(0., 0., 0.),
                          const Vector3D& sdir  = Vector3D(0., 0., 0.),
                          double          spath = 0.)
      : startPos(spos), startDir(sdir), pathLength(spath)
    {
    }

    /// Boolean() operator - returns false for void modifier
    explicit operator bool() const { return true; }

    /// A corrector for step estimation
    ///
    /// @param pos[in,out] the position for the path intersection
    /// @param dir[in,out] the direction for the path intersection
    /// @param path[in,out] path that as a first estimate
    bool
    operator()(Vector3D& pos, Vector3D& dir, double& path) const
    {
      // Approximation with straight line
      if (path * path < straightLineStep * straightLineStep) {
        return false;
      }
      // no correction at start positions
      if (pathLength == 0. || pos.isApprox(startPos)
          || dir.isApprox(startDir)) {
        // can not correct: this is the initial step, no reference point
        return false;
      }
      // change of direction per path length unit
      Vector3D deltaDir = (dir - startDir) * 1. / pathLength;
      dir               = dir + 0.5 * path * deltaDir;
      // return true
      return true;
    }

    /// Step size manipulation call
    bool
    operator()(ConstrainedStep& step) const
    {
      // Don't do anything if you are under accuracy or user control
      auto stepType = step.currentType();
      if (stepType == Cstep::accuracy or stepType == Cstep::user) {
        return false;
      }
      // Apply the step modification
      // Update navigation and target step
      double navStep = stepModification * step.value(Cstep::actor);
      // We need also modify the target stepping for the moment
      step.update(navStep, Cstep::actor, false);
      return true;
    }
  };

}  // namespace detail
}  // namespace Acts
