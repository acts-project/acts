// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_PROPAGATOR_STEP_CALLER_HPP
#define ACTS_PROPAGATOR_STEP_CALLER_HPP 1

#include "ACTS/Extrapolation/Direction.hpp"

namespace Acts {

namespace detail {

  template <typename Stepper, typename Cache, Direction dir>
  struct step_caller;

  template <typename Stepper, typename Cache>
  struct step_caller<Stepper, Cache, forward>
  {
    static double
    step(const Stepper& s, Cache& c, double& stepMax)
    {
      return s.step_forward(c, stepMax);
    }
  };

  template <typename Stepper, typename Cache>
  struct step_caller<Stepper, Cache, backward>
  {
    static double
    step(const Stepper& s, Cache& c, double& stepMax)
    {
      return s.step_backward(c, stepMax);
    }
  };
}  // namespace detail

}  // namespace Acts
#endif  // ACTS_PROPAGATOR_STEP_CALLER_HPP
