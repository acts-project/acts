// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/StraightLineStepper.hpp"

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"

namespace Acts {

auto StraightLineStepper::curvilinearState(State& state, bool transportCov) const
    -> CurvilinearState {
  FreeVector parameters;
  parameters[eFreePos0] = state.pos[ePos0];
  parameters[eFreePos1] = state.pos[ePos1];
  parameters[eFreePos2] = state.pos[ePos2];
  parameters[eFreeTime] = state.t;
  parameters[eFreeDir0] = state.dir[eMom0];
  parameters[eFreeDir1] = state.dir[eMom1];
  parameters[eFreeDir2] = state.dir[eMom2];
  parameters[eFreeQOverP] = (state.q != 0. ? state.q : 1.) / state.p;
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.jacDirToAngle, state.jacAngleToDir, parameters,
      state.covTransport && transportCov, state.pathAccumulated);
}

auto StraightLineStepper::freeState(State& state) const -> FreeState {
  FreeVector parameters;
  parameters[eFreePos0] = state.pos[ePos0];
  parameters[eFreePos1] = state.pos[ePos1];
  parameters[eFreePos2] = state.pos[ePos2];
  parameters[eFreeTime] = state.t;
  parameters[eFreeDir0] = state.dir[eMom0];
  parameters[eFreeDir1] = state.dir[eMom1];
  parameters[eFreeDir2] = state.dir[eMom2];
  parameters[eFreeQOverP] = (state.q != 0. ? state.q : 1.) / state.p;
  return detail::freeState(state.cov, state.jacobian, state.jacTransport,
                           state.derivative, state.jacToGlobal,
                           state.jacDirToAngle, state.jacAngleToDir, parameters,
                           state.covTransport, state.pathAccumulated);
}

auto
StraightLineStepper::boundState(State& state, const Surface& surface,
                                bool transportCov) const -> BoundState {
  FreeVector parameters;
  parameters[eFreePos0] = state.pos[ePos0];
  parameters[eFreePos1] = state.pos[ePos1];
  parameters[eFreePos2] = state.pos[ePos2];
  parameters[eFreeTime] = state.t;
  parameters[eFreeDir0] = state.dir[eMom0];
  parameters[eFreeDir1] = state.dir[eMom1];
  parameters[eFreeDir2] = state.dir[eMom2];
  parameters[eFreeQOverP] = (state.q != 0. ? state.q : 1.) / state.p;
  return detail::boundState(state.geoContext, state.cov, state.jacobian,
                            state.jacTransport, state.derivative,
                            state.jacToGlobal, state.jacDirToAngle,
                            state.jacAngleToDir, parameters, state.covTransport && transportCov,
                            state.pathAccumulated, surface);
}

void StraightLineStepper::update(State& state, const FreeVector& parameters,
                                 const Covariance& covariance) const {
  state.pos = parameters.template segment<3>(eFreePos0);
  state.dir = parameters.template segment<3>(eFreeDir0).normalized();
  state.p = std::abs(1. / parameters[eFreeQOverP]);
  state.t = parameters[eFreeTime];
  state.cov = covariance;
}

void StraightLineStepper::update(State& state, const Vector3D& uposition,
                                 const Vector3D& udirection, double up,
                                 double time) const {
  state.pos = uposition;
  state.dir = udirection.normalized();
  state.p = up;
  state.t = time;
}

void StraightLineStepper::covarianceTransport(State& state,
                                              const Surface& surface) const {
  FreeVector parameters;
  parameters[eFreePos0] = state.pos[ePos0];
  parameters[eFreePos1] = state.pos[ePos1];
  parameters[eFreePos2] = state.pos[ePos2];
  parameters[eFreeTime] = state.t;
  parameters[eFreeDir0] = state.dir[eMom0];
  parameters[eFreeDir1] = state.dir[eMom1];
  parameters[eFreeDir2] = state.dir[eMom2];
  parameters[eFreeQOverP] = state.q / state.p;
  detail::covarianceTransport(state.geoContext, state.cov, state.jacobian,
                              state.jacTransport, state.derivative,
                              state.jacToGlobal, state.jacDirToAngle,
                              state.jacAngleToDir, parameters, surface);
}

void StraightLineStepper::resetState(State& state,
                                     const BoundVector& boundParams,
                                     const BoundSymMatrix& cov,
                                     const Surface& surface,
                                     const NavigationDirection navDir,
                                     const double stepSize) const {
  // Update the stepping state
  update(state,
         detail::transformBoundToFreeParameters(surface, state.geoContext,
                                                boundParams),
         cov);
  state.navDir = navDir;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  surface.initJacobianToGlobal(state.geoContext, state.jacToGlobal,
                               position(state), direction(state), boundParams);
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

}  // namespace Acts
