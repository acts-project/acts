// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename B, typename E, typename A>
Acts::EigenStepper<B, E, A>::EigenStepper(B bField)
    : m_bField(std::move(bField)) {}

template <typename B, typename E, typename A>
auto Acts::EigenStepper<B, E, A>::boundState(State& state,
                                             const Surface& surface,
                                             bool reinitialize) const
    -> BoundState {
  // Transport the covariance to here
  std::optional<Covariance> covOpt = std::nullopt;
  if (state.covTransport) {
    covarianceTransport(state, surface, reinitialize);
    covOpt = state.cov;
  }
  // Create the bound parameters
  BoundParameters parameters(state.geoContext, std::move(covOpt), state.pos,
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

template <typename B, typename E, typename A>
auto Acts::EigenStepper<B, E, A>::curvilinearState(State& state,
                                                   bool reinitialize) const
    -> CurvilinearState {
  // Transport the covariance to here
  std::optional<Covariance> covOpt = std::nullopt;
  if (state.covTransport) {
    covarianceTransport(state, reinitialize);
    covOpt = state.cov;
  }
  // Create the curvilinear parameters
  CurvilinearParameters parameters(std::move(state.cov), state.pos,
                                   state.p * state.dir, state.q, state.t);
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

template <typename B, typename E, typename A>
void Acts::EigenStepper<B, E, A>::update(State& state,
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

template <typename B, typename E, typename A>
void Acts::EigenStepper<B, E, A>::update(State& state,
                                         const Vector3D& uposition,
                                         const Vector3D& udirection, double up,
                                         double time) const {
  state.pos = uposition;
  state.dir = udirection;
  state.p = up;
  state.t = time;
}

template <typename B, typename E, typename A>
void Acts::EigenStepper<B, E, A>::covarianceTransport(State& state,
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
  jacToCurv(5, 3) = 1;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 6) = -invSinTheta;
  jacToCurv(4, 7) = 1;
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
    state.jacToGlobal(3, eT) = 1.;
    state.jacToGlobal(4, ePHI) = -sinTheta * sinPhi;
    state.jacToGlobal(4, eTHETA) = cosTheta * cosPhi;
    state.jacToGlobal(5, ePHI) = sinTheta * cosPhi;
    state.jacToGlobal(5, eTHETA) = cosTheta * sinPhi;
    state.jacToGlobal(6, eTHETA) = -sinTheta;
    state.jacToGlobal(7, eQOP) = 1.;
  }
  // Store The global and bound jacobian (duplication for the moment)
  state.jacobian = jacFull * state.jacobian;
}

template <typename B, typename E, typename A>
void Acts::EigenStepper<B, E, A>::covarianceTransport(State& state,
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

template <typename B, typename E, typename A>
template <typename propagator_state_t>
Acts::Result<double> Acts::EigenStepper<B, E, A>::step(
    propagator_state_t& state) const {
  using namespace UnitLiterals;

  // Runge-Kutta integrator state
  auto& sd = state.stepping.stepData;
  double error_estimate = 0.;
  double h2, half_h;

  // First Runge-Kutta point (at current position)
  sd.B_first = getField(state.stepping, state.stepping.pos);
  if (!state.stepping.extension.validExtensionForStep(state, *this) ||
      !state.stepping.extension.k1(state, *this, sd.k1, sd.B_first, sd.kQoP)) {
    return 0.;
  }

  // The following functor starts to perform a Runge-Kutta step of a certain
  // size, going up to the point where it can return an estimate of the local
  // integration error. The results are stated in the local variables above,
  // allowing integration to continue once the error is deemed satisfactory
  const auto tryRungeKuttaStep = [&](const ConstrainedStep& h) -> bool {
    // State the square and half of the step size
    h2 = h * h;
    half_h = h * 0.5;

    // Second Runge-Kutta point
    const Vector3D pos1 =
        state.stepping.pos + half_h * state.stepping.dir + h2 * 0.125 * sd.k1;
    sd.B_middle = getField(state.stepping, pos1);
    if (!state.stepping.extension.k2(state, *this, sd.k2, sd.B_middle, sd.kQoP,
                                     half_h, sd.k1)) {
      return false;
    }

    // Third Runge-Kutta point
    if (!state.stepping.extension.k3(state, *this, sd.k3, sd.B_middle, sd.kQoP,
                                     half_h, sd.k2)) {
      return false;
    }

    // Last Runge-Kutta point
    const Vector3D pos2 =
        state.stepping.pos + h * state.stepping.dir + h2 * 0.5 * sd.k3;
    sd.B_last = getField(state.stepping, pos2);
    if (!state.stepping.extension.k4(state, *this, sd.k4, sd.B_last, sd.kQoP, h,
                                     sd.k3)) {
      return false;
    }

    // Compute and check the local integration error estimate
    error_estimate = std::max(
        h2 * ((sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>() +
              std::abs(sd.kQoP[0] - sd.kQoP[1] - sd.kQoP[2] + sd.kQoP[3])),
        1e-20);
    return (error_estimate <= state.options.tolerance) &&
           ((h.currentType() != ConstrainedStep::accuracy) ||
            (error_estimate >= state.options.tolerance / 10));
  };

  double stepSizeScaling = 1.;
  size_t nStepTrials = 0;
  // Select and adjust the appropriate Runge-Kutta step size as given
  // ATL-SOFT-PUB-2009-001
  while (!tryRungeKuttaStep(state.stepping.stepSize)) {
    stepSizeScaling =
        std::min(std::max(0.25, std::pow((state.options.tolerance /
                                          std::abs(2. * error_estimate)),
                                         0.25)),
                 4.);
    if (stepSizeScaling == 1.) {
      break;
    }
    state.stepping.stepSize = state.stepping.stepSize * stepSizeScaling;

    // If step size becomes too small the particle remains at the initial
    // place
    if (state.stepping.stepSize * state.stepping.stepSize <
        state.options.stepSizeCutOff * state.options.stepSizeCutOff) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > state.options.maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
    nStepTrials++;
  }

  // use the adjusted step size
  const double h = state.stepping.stepSize;

  // When doing error propagation, update the associated Jacobian matrix
  if (state.stepping.covTransport) {
    // The step transport matrix in global coordinates
    FreeMatrix D;
    if (!state.stepping.extension.finalize(state, *this, h, D)) {
      return EigenStepperError::StepInvalid;
    }

    // for moment, only update the transport part
    state.stepping.jacTransport = D * state.stepping.jacTransport;
  } else {
    if (!state.stepping.extension.finalize(state, *this, h)) {
      return EigenStepperError::StepInvalid;
    }
  }

  // Update the track parameters according to the equations of motion
  state.stepping.pos +=
      h * state.stepping.dir + h2 / 6. * (sd.k1 + sd.k2 + sd.k3);
  state.stepping.dir += h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
  state.stepping.dir /= state.stepping.dir.norm();
  if (state.stepping.covTransport) {
    state.stepping.derivative.template head<3>() = state.stepping.dir;
    state.stepping.derivative.template segment<3>(4) = sd.k4;
  }
  state.stepping.pathAccumulated += h;
  return h;
}