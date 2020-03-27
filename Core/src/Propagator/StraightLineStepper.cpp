// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/StraightLineStepper.hpp"

namespace Acts {

std::tuple<BoundParameters, BoundMatrix, double>
StraightLineStepper::boundState(State& state, const Surface& surface,
                                bool reinitialize) const {
  // Transport the covariance to here
  std::optional<Covariance> cov = std::nullopt;
  if (state.covTransport) {
    covarianceTransport(state, surface, reinitialize);
    cov = state.cov;
  }
  // Create the bound parameters
  BoundParameters parameters(state.geoContext, cov, state.pos,
                             state.p * state.dir, state.q, state.t,
                             surface.getSharedPtr());
  // Create the bound state
  BoundState bState{std::move(parameters), state.jacobian,
                    state.pathAccumulated};
  // Reset the jacobian to identity
  if (reinitialize) {
    state.jacobian = Jacobian::Identity();
  }
  /// Return the State
  return bState;
}

std::tuple<CurvilinearParameters, BoundMatrix, double>
StraightLineStepper::curvilinearState(State& state, bool reinitialize) const {
  // Transport the covariance to here
  std::optional<Covariance> cov = std::nullopt;
  if (state.covTransport) {
    covarianceTransport(state, reinitialize);
    cov = state.cov;
  }
  // Create the curvilinear parameters
  CurvilinearParameters parameters(cov, state.pos, state.p * state.dir, state.q,
                                   state.t);
  // Create the bound state
  CurvilinearState curvState{std::move(parameters), state.jacobian,
                             state.pathAccumulated};
  // Reset the jacobian to identity
  if (reinitialize) {
    state.jacobian = Jacobian::Identity();
  }
  /// Return the State
  return curvState;
}

void StraightLineStepper::update(State& state,
                                 const BoundParameters& pars) const {
  const auto& mom = pars.momentum();
  state.pos = pars.position();
  state.dir = mom.normalized();
  state.p = mom.norm();
  state.t = pars.time();

  if (pars.covariance()) {
    state.cov = (*(pars.covariance()));
  }
}

void StraightLineStepper::update(State& state, const Vector3D& uposition,
                                 const Vector3D& udirection, double up,
                                 double time) const {
  state.pos = uposition;
  state.dir = udirection;
  state.p = up;
  state.t = time;
}

void StraightLineStepper::covarianceTransport(State& state,
                                              bool reinitialize) const {
  // Optimized trigonometry on the propagation direction
  const double x = state.dir(0);  // == cos(phi) * sin(theta)
  const double y = state.dir(1);  // == sin(phi) * sin(theta)
  const double z = state.dir(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // prepare the jacobian to curvilinear
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(0, 0) = -sinPhi;
    jacToCurv(0, 1) = cosPhi;
    jacToCurv(1, 0) = -cosPhi * cosTheta;
    jacToCurv(1, 1) = -sinPhi * cosTheta;
    jacToCurv(1, 2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = sqrt(y * y + z * z);
    const double invC = 1. / c;
    jacToCurv(0, 1) = -z * invC;
    jacToCurv(0, 2) = y * invC;
    jacToCurv(1, 0) = c;
    jacToCurv(1, 1) = -x * y * invC;
    jacToCurv(1, 2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(5, 3) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 6) = -invSinTheta;
  jacToCurv(4, 7) = 1.;
  // Apply the transport from the steps on the jacobian
  state.jacToGlobal = state.jacTransport * state.jacToGlobal;
  // Transport the covariance
  ActsRowVectorD<3> normVec(state.dir);
  const BoundRowVector sfactors =
      normVec *
      state.jacToGlobal.template topLeftCorner<3, eBoundParametersSize>();
  // The full jacobian is ([to local] jacobian) * ([transport] jacobian)
  const Jacobian jacFull =
      jacToCurv * (state.jacToGlobal - state.derivative * sfactors);
  // Apply the actual covariance transport
  state.cov = (jacFull * state.cov * jacFull.transpose());
  // Reinitialize if asked to do so
  // this is useful for interruption calls
  if (reinitialize) {
    // reset the jacobians
    state.jacToGlobal = BoundToFreeMatrix::Zero();
    state.jacTransport = FreeMatrix::Identity();
    // fill the jacobian to global for next transport
    state.jacToGlobal(0, eLOC_0) = -sinPhi;
    state.jacToGlobal(0, eLOC_1) = -cosPhi * cosTheta;
    state.jacToGlobal(1, eLOC_0) = cosPhi;
    state.jacToGlobal(1, eLOC_1) = -sinPhi * cosTheta;
    state.jacToGlobal(2, eLOC_1) = sinTheta;
    state.jacToGlobal(3, eT) = 1;
    state.jacToGlobal(4, ePHI) = -sinTheta * sinPhi;
    state.jacToGlobal(4, eTHETA) = cosTheta * cosPhi;
    state.jacToGlobal(5, ePHI) = sinTheta * cosPhi;
    state.jacToGlobal(5, eTHETA) = cosTheta * sinPhi;
    state.jacToGlobal(6, eTHETA) = -sinTheta;
    state.jacToGlobal(7, eQOP) = 1;
  }
  // Store The global and bound jacobian (duplication for the moment)
  state.jacobian = jacFull * state.jacobian;
}

void StraightLineStepper::covarianceTransport(State& state,
                                              const Surface& surface,
                                              bool reinitialize) const {
  using VectorHelpers::phi;
  using VectorHelpers::theta;
  // Initialize the transport final frame jacobian
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // initalize the jacobian to local, returns the transposed ref frame
  auto rframeT = surface.initJacobianToLocal(state.geoContext, jacToLocal,
                                             state.pos, state.dir);
  // Update the jacobian with the transport from the steps
  state.jacToGlobal = state.jacTransport * state.jacToGlobal;
  // calculate the form factors for the derivatives
  const BoundRowVector sVec = surface.derivativeFactors(
      state.geoContext, state.pos, state.dir, rframeT, state.jacToGlobal);
  // the full jacobian is ([to local] jacobian) * ([transport] jacobian)
  const Jacobian jacFull =
      jacToLocal * (state.jacToGlobal - state.derivative * sVec);
  // Apply the actual covariance transport
  state.cov = (jacFull * state.cov * jacFull.transpose());
  // Reinitialize if asked to do so
  // this is useful for interruption calls
  if (reinitialize) {
    // reset the jacobians
    state.jacToGlobal = BoundToFreeMatrix::Zero();
    state.jacTransport = FreeMatrix::Identity();
    // reset the derivative
    state.derivative = FreeVector::Zero();
    // fill the jacobian to global for next transport
    Vector2D loc{0., 0.};
    surface.globalToLocal(state.geoContext, state.pos, state.dir, loc);
    BoundVector pars;
    pars << loc[eLOC_0], loc[eLOC_1], phi(state.dir), theta(state.dir),
        state.q / state.p, state.t;
    surface.initJacobianToGlobal(state.geoContext, state.jacToGlobal, state.pos,
                                 state.dir, pars);
  }
  // Store The global and bound jacobian (duplication for the moment)
  state.jacobian = jacFull * state.jacobian;
}
}  // namespace Acts
