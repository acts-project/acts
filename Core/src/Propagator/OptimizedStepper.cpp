// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/OptimizedStepper.hpp"

#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"

#include <cmath>
#include <cstdint>

namespace Acts {

namespace {

template <typename T>
void rk4_k1(const T* d, const T lambda, const T* B, T* k1) {
  k1[0] = lambda * (-B[1] * d[2] + B[2] * d[1]);
  k1[1] = lambda * (B[0] * d[2] - B[2] * d[0]);
  k1[2] = lambda * (-B[0] * d[1] + B[1] * d[0]);
}

template <typename T>
void rk4_p2(const T* p, const T* d, const T h, const T* k1, T* p2) {
  p2[0] = (1.0 / 8.0) * std::pow(h, 2) * k1[0] + (1.0 / 2.0) * h * d[0] + p[0];
  p2[1] = (1.0 / 8.0) * std::pow(h, 2) * k1[1] + (1.0 / 2.0) * h * d[1] + p[1];
  p2[2] = (1.0 / 8.0) * std::pow(h, 2) * k1[2] + (1.0 / 2.0) * h * d[2] + p[2];
}

template <typename T>
void rk4_k2(const T* d, const T h, const T lambda, const T* B, const T* k1,
            T* k2) {
  const auto x0 = (1.0 / 2.0) * h;
  const auto x1 = x0 * k1[1] + d[1];
  const auto x2 = x0 * k1[2] + d[2];
  const auto x3 = x0 * k1[0] + d[0];
  k2[0] = lambda * (x1 * B[2] - x2 * B[1]);
  k2[1] = lambda * (x2 * B[0] - x3 * B[2]);
  k2[2] = lambda * (-x1 * B[0] + x3 * B[1]);
}

template <typename T>
void rk4_k3(const T* d, const T h, const T lambda, const T* B, const T* k2,
            T* k3) {
  const auto x0 = (1.0 / 2.0) * h;
  const auto x1 = x0 * k2[1] + d[1];
  const auto x2 = x0 * k2[2] + d[2];
  const auto x3 = x0 * k2[0] + d[0];
  k3[0] = lambda * (x1 * B[2] - x2 * B[1]);
  k3[1] = lambda * (x2 * B[0] - x3 * B[2]);
  k3[2] = lambda * (-x1 * B[0] + x3 * B[1]);
}

template <typename T>
void rk4_p3(const T* p, const T* d, const T h, const T* k3, T* p3) {
  p3[0] = (1.0 / 2.0) * std::pow(h, 2) * k3[0] + h * d[0] + p[0];
  p3[1] = (1.0 / 2.0) * std::pow(h, 2) * k3[1] + h * d[1] + p[1];
  p3[2] = (1.0 / 2.0) * std::pow(h, 2) * k3[2] + h * d[2] + p[2];
}

template <typename T>
void rk4_k4(const T* d, const T h, const T lambda, const T* B, const T* k3,
            T* k4) {
  const auto x0 = h * k3[1] + d[1];
  const auto x1 = h * k3[2] + d[2];
  const auto x2 = h * k3[0] + d[0];
  k4[0] = lambda * (x0 * B[2] - x1 * B[1]);
  k4[1] = lambda * (x1 * B[0] - x2 * B[2]);
  k4[2] = lambda * (-x0 * B[0] + x2 * B[1]);
}

template <typename T>
void rk4_err(const T h, const T* k1, const T* k2, const T* k3, const T* k4,
             T* err) {
  *err = std::pow(h, 2) * (std::abs(k1[0] - k2[0] - k3[0] + k4[0]) +
                           std::abs(k1[1] - k2[1] - k3[1] + k4[1]) +
                           std::abs(k1[2] - k2[2] - k3[2] + k4[2]));
}

template <typename T>
void rk4_fin(const T* p, const T* d, const T h, const T* k1, const T* k2,
             const T* k3, const T* k4, T* new_p, T* new_d) {
  const auto x0 = (1.0 / 6.0) * h;
  new_p[0] =
      (1.0 / 6.0) * std::pow(h, 2) * (k1[0] + k2[0] + k3[0]) + h * d[0] + p[0];
  new_p[1] =
      (1.0 / 6.0) * std::pow(h, 2) * (k1[1] + k2[1] + k3[1]) + h * d[1] + p[1];
  new_p[2] =
      (1.0 / 6.0) * std::pow(h, 2) * (k1[2] + k2[2] + k3[2]) + h * d[2] + p[2];
  new_d[0] =
      (x0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) + d[0]) /
      std::sqrt(
          std::pow(
              std::abs(x0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) + d[0]),
              2) +
          std::pow(
              std::abs(x0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) + d[1]),
              2) +
          std::pow(
              std::abs(x0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) + d[2]),
              2));
  new_d[1] =
      (x0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) + d[1]) /
      std::sqrt(
          std::pow(
              std::abs(x0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) + d[0]),
              2) +
          std::pow(
              std::abs(x0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) + d[1]),
              2) +
          std::pow(
              std::abs(x0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) + d[2]),
              2));
  new_d[2] =
      (x0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) + d[2]) /
      std::sqrt(
          std::pow(
              std::abs(x0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) + d[0]),
              2) +
          std::pow(
              std::abs(x0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) + d[1]),
              2) +
          std::pow(
              std::abs(x0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) + d[2]),
              2));
}

}  // namespace

OptimizedStepper::OptimizedStepper(
    std::shared_ptr<const MagneticFieldProvider> bField, double overstepLimit)
    : m_bField(std::move(bField)), m_overstepLimit(overstepLimit) {}

OptimizedStepper::State OptimizedStepper::makeState(
    std::reference_wrapper<const GeometryContext> gctx,
    std::reference_wrapper<const MagneticFieldContext> mctx,
    const BoundTrackParameters& par, double ssize) const {
  return State{gctx, m_bField->makeCache(mctx), par, ssize};
}

void OptimizedStepper::resetState(State& state, const BoundVector& boundParams,
                                  const BoundSquareMatrix& cov,
                                  const Surface& surface,
                                  const double stepSize) const {
  FreeVector freeParams =
      transformBoundToFreeParameters(surface, state.geoContext, boundParams);

  // Update the stepping state
  state.pars = freeParams;
  state.cov = cov;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

Result<std::tuple<BoundTrackParameters, BoundMatrix, double>>
OptimizedStepper::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  return detail::boundState(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated,
      freeToBoundCorrection);
}

std::tuple<CurvilinearTrackParameters, BoundMatrix, double>
OptimizedStepper::curvilinearState(State& state, bool transportCov) const {
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated);
}

void OptimizedStepper::update(State& state, const FreeVector& freeParams,
                              const BoundVector& /*boundParams*/,
                              const Covariance& covariance,
                              const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
}

void OptimizedStepper::update(State& state, const Vector3& uposition,
                              const Vector3& udirection, double qop,
                              double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qop;
}

void OptimizedStepper::transportCovarianceToCurvilinear(State& state) const {
  detail::transportCovarianceToCurvilinear(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars.template segment<3>(eFreeDir0));
}

void OptimizedStepper::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::transportCovarianceToBound(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, freeToBoundCorrection);
}

Result<double> OptimizedStepper::stepImpl(
    State& state, Direction stepDirection, double stepTolerance,
    double stepSizeCutOff, std::size_t maxRungeKuttaStepTrials) const {
  // Runge-Kutta integrator state
  auto& sd = state.stepData;

  auto pos = position(state);
  auto dir = direction(state);
  auto qop = qOverP(state);
  auto m = particleHypothesis(state).mass();
  auto p = absoluteMomentum(state);

  // First Runge-Kutta point (at current position)
  auto fieldRes = getField(state, pos);
  if (!fieldRes.ok()) {
    return fieldRes.error();
  }
  sd.B_first = *fieldRes;
  rk4_k1(dir.data(), qop, sd.B_first.data(), sd.k1.data());

  double h = state.stepSize.value() * stepDirection;
  std::size_t nStepTrials = 0;
  double error_estimate = 0.;

  while (true) {
    // Second Runge-Kutta point
    Vector3 pos2;
    rk4_p2(pos.data(), dir.data(), h, sd.k1.data(), pos2.data());
    auto field = getField(state, pos2);
    if (!field.ok()) {
      return field.error();
    }
    sd.B_middle = *field;
    rk4_k2(dir.data(), h, qop, sd.B_middle.data(), sd.k1.data(), sd.k2.data());

    // Third Runge-Kutta point
    rk4_k3(dir.data(), h, qop, sd.B_middle.data(), sd.k2.data(), sd.k3.data());

    // Last Runge-Kutta point
    Vector3 pos3;
    rk4_p3(pos.data(), dir.data(), h, sd.k1.data(), pos2.data());
    field = getField(state, pos3);
    if (!field.ok()) {
      return field.error();
    }
    sd.B_last = *field;
    rk4_k4(dir.data(), h, qop, sd.B_last.data(), sd.k3.data(), sd.k4.data());

    // Compute and check the local integration error estimate
    rk4_err(h, sd.k1.data(), sd.k2.data(), sd.k3.data(), sd.k4.data(),
            &error_estimate);

    if (error_estimate <= stepTolerance) {
      break;
    }

    h *= 0.5;
    state.stepSize.setAccuracy(h);

    // If step size becomes too small the particle remains at the initial
    // place
    if (std::abs(h) < std::abs(stepSizeCutOff)) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
    nStepTrials++;
  }

  state.pathAccumulated += h;
  state.stepSize.nStepTrials = nStepTrials;

  // Update the track parameters according to the equations of motion
  rk4_fin(pos.data(), dir.data(), h, sd.k1.data(), sd.k2.data(), sd.k3.data(),
          sd.k4.data(), state.pars.template segment<3>(eFreeDir0).data(),
          state.pars.template segment<3>(eFreeDir0).data());

  /// This evaluation is based on dt/ds = 1/v = 1/(beta * c) with the velocity
  /// v, the speed of light c and beta = v/c. This can be re-written as dt/ds
  /// = sqrt(m^2/p^2 + c^{-2}) with the mass m and the momentum p.
  auto dtds = std::sqrt(1 + m * m / (p * p));
  state.pars[eFreeTime] += h * dtds;

  if (state.covTransport) {
    double half_h = 0.5 * h;

    // When doing error propagation, update the associated Jacobian matrix
    state.derivative.template head<3>() =
        state.pars.template segment<3>(eFreeDir0);
    state.derivative[3] = dtds;
    state.derivative.template segment<3>(4) = sd.k4;

    // The step transport matrix in global coordinates
    FreeMatrix D = FreeMatrix::Identity();

    // This sets the reference to the sub matrices
    // dFdx is already initialised as (3x3) idendity
    auto dFdT = D.block<3, 3>(0, 4);
    auto dFdL = D.block<3, 1>(0, 7);
    // dGdx is already initialised as (3x3) zero
    auto dGdT = D.block<3, 3>(4, 4);
    auto dGdL = D.block<3, 1>(4, 7);

    ActsMatrix<3, 3> dk1dT = ActsMatrix<3, 3>::Zero();
    ActsMatrix<3, 3> dk2dT = ActsMatrix<3, 3>::Identity();
    ActsMatrix<3, 3> dk3dT = ActsMatrix<3, 3>::Identity();
    ActsMatrix<3, 3> dk4dT = ActsMatrix<3, 3>::Identity();

    Vector3 dk1dL = Vector3::Zero();
    Vector3 dk2dL = Vector3::Zero();
    Vector3 dk3dL = Vector3::Zero();
    Vector3 dk4dL = Vector3::Zero();

    // For the case without energy loss
    dk1dL = dir.cross(sd.B_first);
    dk2dL = (dir + half_h * sd.k1).cross(sd.B_middle) +
            qop * half_h * dk1dL.cross(sd.B_middle);
    dk3dL = (dir + half_h * sd.k2).cross(sd.B_middle) +
            qop * half_h * dk2dL.cross(sd.B_middle);
    dk4dL =
        (dir + h * sd.k3).cross(sd.B_last) + qop * h * dk3dL.cross(sd.B_last);

    dk1dT(0, 1) = sd.B_first.z();
    dk1dT(0, 2) = -sd.B_first.y();
    dk1dT(1, 0) = -sd.B_first.z();
    dk1dT(1, 2) = sd.B_first.x();
    dk1dT(2, 0) = sd.B_first.y();
    dk1dT(2, 1) = -sd.B_first.x();
    dk1dT *= qop;

    dk2dT += half_h * dk1dT;
    dk2dT = qop * VectorHelpers::cross(dk2dT, sd.B_middle);

    dk3dT += half_h * dk2dT;
    dk3dT = qop * VectorHelpers::cross(dk3dT, sd.B_middle);

    dk4dT += h * dk3dT;
    dk4dT = qop * VectorHelpers::cross(dk4dT, sd.B_last);

    dFdT.setIdentity();
    dFdT += h / 6. * (dk1dT + dk2dT + dk3dT);
    dFdT *= h;

    dFdL = (h * h) / 6. * (dk1dL + dk2dL + dk3dL);

    dGdT += h / 6. * (dk1dT + 2. * (dk2dT + dk3dT) + dk4dT);

    dGdL = h / 6. * (dk1dL + 2. * (dk2dL + dk3dL) + dk4dL);

    D(3, 7) = h * m * m * qop / dtds;

    // See the documentation of Acts::blockedMult for a description of blocked
    // matrix multiplication. However, we can go one further. Let's assume that
    // some of these sub-matrices are zero matrices 0₈ and identity matrices
    // I₈, namely:
    //
    // D₁₁ = I₈, J₁₁ = I₈, D₂₁ = 0₈, J₂₁ = 0₈
    //
    // Which gives:
    //
    // K₁₁ = I₈  * I₈  + D₁₂ * 0₈  = I₈
    // K₁₂ = I₈  * J₁₂ + D₁₂ * J₂₂ = J₁₂ + D₁₂ * J₂₂
    // K₂₁ = 0₈  * I₈  + D₂₂ * 0₈  = 0₈
    // K₂₂ = 0₈  * J₁₂ + D₂₂ * J₂₂ = D₂₂ * J₂₂
    //
    // Furthermore, we're constructing K in place of J, and since
    // K₁₁ = I₈ = J₁₁ and K₂₁ = 0₈ = D₂₁, we don't actually need to touch those
    // sub-matrices at all!
    assert((D.topLeftCorner<4, 4>().isIdentity()));
    assert((D.bottomLeftCorner<4, 4>().isZero()));
    assert((state.jacTransport.template topLeftCorner<4, 4>().isIdentity()));
    assert((state.jacTransport.template bottomLeftCorner<4, 4>().isZero()));

    state.jacTransport.template topRightCorner<4, 4>() +=
        D.topRightCorner<4, 4>() *
        state.jacTransport.template bottomRightCorner<4, 4>();
    state.jacTransport.template bottomRightCorner<4, 4>() =
        (D.bottomRightCorner<4, 4>() *
         state.jacTransport.template bottomRightCorner<4, 4>())
            .eval();
  }

  return h;
}

void OptimizedStepper::setIdentityJacobian(State& state) const {
  state.jacobian = BoundMatrix::Identity();
}

}  // namespace Acts
