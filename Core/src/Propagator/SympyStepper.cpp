// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/SympyStepper.hpp"

#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"

#include <cmath>
#include <cstdint>

namespace Acts {

namespace {

template <typename T, typename GetB>
bool rk4(const T* p, const T* d, const T h, const T lambda, const T m,
         const T p_abs, T* new_p, T* new_d, T* new_time, T* J, GetB getB,
         bool covTransport) {
  const auto x7 = std::pow(h, 2);
  const auto x11 = (1.0 / 2.0) * h;
  const auto B1 = getB(p);
  const auto x6 = h * d[0];
  const auto x9 = h * d[1];
  const auto x10 = h * d[2];
  const auto x8 = (1.0 / 8.0) * x7;
  const auto x24 = (1.0 / 2.0) * x7;
  const auto x0 = -B1[1] * d[2] + B1[2] * d[1];
  const auto x2 = B1[0] * d[2] - B1[2] * d[0];
  const auto x4 = -B1[0] * d[1] + B1[1] * d[0];
  const auto x25 = x6 + p[0];
  const auto x26 = x9 + p[1];
  const auto x27 = x10 + p[2];
  const auto x1 = lambda * x0;
  const auto x3 = lambda * x2;
  const auto x5 = lambda * x4;
  T k1[3];
  k1[0] = x1;
  k1[1] = x3;
  k1[2] = x5;
  const auto x15 = x11 * k1[0] + d[0];
  const auto x12 = x11 * k1[1] + d[1];
  const auto x13 = x11 * k1[2] + d[2];
  T p2[3];
  p2[0] = (1.0 / 2.0) * x6 + x8 * k1[0] + p[0];
  p2[1] = (1.0 / 2.0) * x9 + x8 * k1[1] + p[1];
  p2[2] = (1.0 / 2.0) * x10 + x8 * k1[2] + p[2];
  const auto B2 = getB(p2);
  const auto x17 = x12 * B2[0];
  const auto x14 = x13 * B2[1];
  const auto x16 = x15 * B2[2];
  T k2[3];
  k2[0] = lambda * (x12 * B2[2] - x14);
  k2[1] = lambda * (x13 * B2[0] - x16);
  k2[2] = lambda * (x15 * B2[1] - x17);
  const auto x21 = x11 * k2[0] + d[0];
  const auto x18 = x11 * k2[1] + d[1];
  const auto x19 = x11 * k2[2] + d[2];
  const auto x23 = x18 * B2[0];
  const auto x20 = x19 * B2[1];
  const auto x22 = x21 * B2[2];
  T k3[3];
  k3[0] = lambda * (x18 * B2[2] - x20);
  k3[1] = lambda * (x19 * B2[0] - x22);
  k3[2] = lambda * (x21 * B2[1] - x23);
  const auto x31 = h * k3[0] + d[0];
  const auto x28 = h * k3[1] + d[1];
  const auto x29 = h * k3[2] + d[2];
  T p3[3];
  p3[0] = x24 * k3[0] + x25;
  p3[1] = x24 * k3[1] + x26;
  p3[2] = x24 * k3[2] + x27;
  const auto B3 = getB(p3);
  const auto x33 = x28 * B3[0];
  const auto x30 = x29 * B3[1];
  const auto x32 = x31 * B3[2];
  T k4[3];
  k4[0] = lambda * (x28 * B3[2] - x30);
  k4[1] = lambda * (x29 * B3[0] - x32);
  k4[2] = lambda * (x31 * B3[1] - x33);
  const auto x34 = k1[0] + k4[0];
  const auto x35 = k1[1] + k4[1];
  const auto x36 = k1[2] + k4[2];
  const auto err =
      x7 * (std::fabs(-x34 + k2[0] + k3[0]) + std::fabs(-x35 + k2[1] + k3[1]) +
            std::fabs(-x36 + k2[2] + k3[2]));
  if (err > 1e-4) {
    return false;
  }
  const auto x37 = (1.0 / 6.0) * x7;
  new_p[0] = x25 + x37 * (k1[0] + k2[0] + k3[0]);
  new_p[1] = x26 + x37 * (k1[1] + k2[1] + k3[1]);
  new_p[2] = x27 + x37 * (k1[2] + k2[2] + k3[2]);
  const auto x38 = (1.0 / 6.0) * h;
  new_d[0] = x38 * (x34 + 2 * k2[0] + 2 * k3[0]) + d[0];
  new_d[1] = x38 * (x35 + 2 * k2[1] + 2 * k3[1]) + d[1];
  new_d[2] = x38 * (x36 + 2 * k2[2] + 2 * k3[2]) + d[2];

  /// This evaluation is based on dt/ds = 1/v = 1/(beta * c) with the velocity
  /// v, the speed of light c and beta = v/c. This can be re-written as dt/ds
  /// = sqrt(m^2/p^2 + c^{-2}) with the mass m and the momentum p.
  auto dtds = std::sqrt(1 + m * m / (p_abs * p_abs));
  *new_time += h * dtds;

  if (covTransport) {
    const auto x75 = (1.0 / 3.0) * h;
    const auto x39 = std::pow(lambda, 2) * x11;
    const auto x48 = lambda * B2[0];
    const auto x45 = lambda * B2[1];
    const auto x44 = lambda * B2[2];
    const auto x56 = lambda * B3[0];
    const auto x52 = lambda * B3[1];
    const auto x54 = lambda * B3[2];
    const auto x58 = x11 * x5;
    const auto x60 = lambda * x37;
    const auto x77 = lambda * x38;
    const auto x46 = x39 * B2[0];
    const auto x40 = x39 * B2[1];
    const auto x42 = x39 * B2[2];
    const auto x49 = x11 * x45;
    const auto x50 = x11 * x44;
    const auto x51 = x11 * x48;
    const auto x53 = h * x52;
    const auto x55 = h * x54;
    const auto x57 = h * x56;
    const auto x67 = x60 * B1[0];
    const auto x63 = x60 * B1[1];
    const auto x61 = x60 * B1[2];
    const auto x84 = x77 * B1[0];
    const auto x80 = x77 * B1[1];
    const auto x78 = x77 * B1[2];
    T dk2dL[3];
    dk2dL[0] =
        (1.0 / 2.0) * h * lambda * x2 * B2[2] + x12 * B2[2] - x14 - x58 * B2[1];
    dk2dL[1] = -x1 * x11 * B2[2] + x13 * B2[0] - x16 + x58 * B2[0];
    dk2dL[2] = (1.0 / 2.0) * h * lambda * x0 * B2[1] - x11 * x3 * B2[0] +
               x15 * B2[1] - x17;
    const auto x47 = x46 * B1[0];
    const auto x41 = x40 * B1[1];
    const auto x43 = x42 * B1[2];
    T dk3dL[3];
    dk3dL[0] = (1.0 / 2.0) * h * lambda * B2[2] * dk2dL[1] + x18 * B2[2] - x20 -
               x49 * dk2dL[2];
    dk3dL[1] = x19 * B2[0] - x22 - x50 * dk2dL[0] + x51 * dk2dL[2];
    dk3dL[2] = (1.0 / 2.0) * h * lambda * B2[1] * dk2dL[0] + x21 * B2[1] - x23 -
               x51 * dk2dL[1];
    const auto x72 = x0 * x37 + x37 * dk2dL[0] + x37 * dk3dL[0];
    const auto x73 = x2 * x37 + x37 * dk2dL[1] + x37 * dk3dL[1];
    const auto x74 = x37 * x4 + x37 * dk2dL[2] + x37 * dk3dL[2];
    T dk2dT[9];
    dk2dT[0] = -x41 - x43;
    dk2dT[1] = x40 * B1[0] + x44;
    dk2dT[2] = x42 * B1[0] - x45;
    dk2dT[3] = -x44 + x46 * B1[1];
    dk2dT[4] = -x43 - x47;
    dk2dT[5] = x42 * B1[1] + x48;
    dk2dT[6] = x45 + x46 * B1[2];
    dk2dT[7] = x40 * B1[2] - x48;
    dk2dT[8] = -x41 - x47;
    T dk4dL[3];
    dk4dL[0] =
        h * lambda * B3[2] * dk3dL[1] + x28 * B3[2] - x30 - x53 * dk3dL[2];
    dk4dL[1] = x29 * B3[0] - x32 - x55 * dk3dL[0] + x57 * dk3dL[2];
    dk4dL[2] =
        h * lambda * B3[1] * dk3dL[0] + x31 * B3[1] - x33 - x57 * dk3dL[1];
    const auto x89 =
        x0 * x38 + x38 * dk4dL[0] + x75 * dk2dL[0] + x75 * dk3dL[0];
    const auto x90 =
        x2 * x38 + x38 * dk4dL[1] + x75 * dk2dL[1] + x75 * dk3dL[1];
    const auto x91 =
        x38 * x4 + x38 * dk4dL[2] + x75 * dk2dL[2] + x75 * dk3dL[2];
    T dk3dT[9];
    dk3dT[0] = (1.0 / 2.0) * h * lambda * B2[2] * dk2dT[3] - x49 * dk2dT[6];
    dk3dT[1] = x44 - x49 * dk2dT[7] + x50 * dk2dT[4];
    dk3dT[2] =
        (1.0 / 2.0) * h * lambda * B2[2] * dk2dT[5] - x45 - x49 * dk2dT[8];
    dk3dT[3] =
        (1.0 / 2.0) * h * lambda * B2[0] * dk2dT[6] - x44 - x50 * dk2dT[0];
    dk3dT[4] = -x50 * dk2dT[1] + x51 * dk2dT[7];
    dk3dT[5] = x48 - x50 * dk2dT[2] + x51 * dk2dT[8];
    dk3dT[6] = x45 + x49 * dk2dT[0] - x51 * dk2dT[3];
    dk3dT[7] =
        (1.0 / 2.0) * h * lambda * B2[1] * dk2dT[1] - x48 - x51 * dk2dT[4];
    dk3dT[8] = (1.0 / 2.0) * h * lambda * B2[1] * dk2dT[2] - x51 * dk2dT[5];
    const auto x59 = x37 * dk2dT[0] + x37 * dk3dT[0] + 1;
    const auto x62 = x37 * dk2dT[1] + x37 * dk3dT[1] + x61;
    const auto x64 = x37 * dk2dT[2] + x37 * dk3dT[2] - x63;
    const auto x65 = x37 * dk2dT[3] + x37 * dk3dT[3] - x61;
    const auto x66 = x37 * dk2dT[4] + x37 * dk3dT[4] + 1;
    const auto x68 = x37 * dk2dT[5] + x37 * dk3dT[5] + x67;
    const auto x69 = x37 * dk2dT[6] + x37 * dk3dT[6] + x63;
    const auto x70 = x37 * dk2dT[7] + x37 * dk3dT[7] - x67;
    const auto x71 = x37 * dk2dT[8] + x37 * dk3dT[8] + 1;
    T dk4dT[9];
    dk4dT[0] = h * lambda * B3[2] * dk3dT[3] - x53 * dk3dT[6];
    dk4dT[1] = -x53 * dk3dT[7] + x54 + x55 * dk3dT[4];
    dk4dT[2] = h * lambda * B3[2] * dk3dT[5] - x52 - x53 * dk3dT[8];
    dk4dT[3] = h * lambda * B3[0] * dk3dT[6] - x54 - x55 * dk3dT[0];
    dk4dT[4] = -x55 * dk3dT[1] + x57 * dk3dT[7];
    dk4dT[5] = -x55 * dk3dT[2] + x56 + x57 * dk3dT[8];
    dk4dT[6] = x52 + x53 * dk3dT[0] - x57 * dk3dT[3];
    dk4dT[7] = h * lambda * B3[1] * dk3dT[1] - x56 - x57 * dk3dT[4];
    dk4dT[8] = h * lambda * B3[1] * dk3dT[2] - x57 * dk3dT[5];
    const auto x76 = x38 * dk4dT[0] + x75 * dk2dT[0] + x75 * dk3dT[0];
    const auto x79 = x38 * dk4dT[1] + x75 * dk2dT[1] + x75 * dk3dT[1] + x78;
    const auto x81 = x38 * dk4dT[2] + x75 * dk2dT[2] + x75 * dk3dT[2] - x80;
    const auto x82 = x38 * dk4dT[3] + x75 * dk2dT[3] + x75 * dk3dT[3] - x78;
    const auto x83 = x38 * dk4dT[4] + x75 * dk2dT[4] + x75 * dk3dT[4];
    const auto x85 = x38 * dk4dT[5] + x75 * dk2dT[5] + x75 * dk3dT[5] + x84;
    const auto x86 = x38 * dk4dT[6] + x75 * dk2dT[6] + x75 * dk3dT[6] + x80;
    const auto x87 = x38 * dk4dT[7] + x75 * dk2dT[7] + x75 * dk3dT[7] - x84;
    const auto x88 = x38 * dk4dT[8] + x75 * dk2dT[8] + x75 * dk3dT[8];

    double j4 = x59 + x76 * J[4] + x82 * J[5] + x86 * J[6];
    double j5 = x62 + x79 * J[4] + x83 * J[5] + x87 * J[6];
    double j6 = x64 + x81 * J[4] + x85 * J[5] + x88 * J[6];
    double j7 = x72 + x89 * J[4] + x90 * J[5] + x91 * J[6];

    double j12 = x65 + x76 * J[12] + x82 * J[13] + x86 * J[14];
    double j13 = x66 + x79 * J[12] + x83 * J[13] + x87 * J[14];
    double j14 = x68 + x81 * J[12] + x85 * J[13] + x88 * J[14];
    double j15 = x73 + x89 * J[12] + x90 * J[13] + x91 * J[14];

    double j20 = x69 + x76 * J[20] + x82 * J[21] + x86 * J[22];
    double j21 = x70 + x79 * J[20] + x83 * J[21] + x87 * J[22];
    double j22 = x71 + x81 * J[20] + x85 * J[21] + x88 * J[22];
    double j23 = x74 + x89 * J[20] + x90 * J[21] + x91 * J[22];

    double j36 = x76 * J[36] + x82 * J[37] + x86 * J[38];
    double j37 = x79 * J[36] + x83 * J[37] + x87 * J[38];
    double j38 = x81 * J[36] + x85 * J[37] + x88 * J[38];
    double j39 = x89 * J[36] + x90 * J[37] + x91 * J[38];

    double j44 = x76 * J[44] + x82 * J[45] + x86 * J[46];
    double j45 = x79 * J[44] + x83 * J[45] + x87 * J[46];
    double j46 = x81 * J[44] + x85 * J[45] + x88 * J[46];
    double j47 = x89 * J[44] + x90 * J[45] + x91 * J[46];

    double j52 = x76 * J[52] + x82 * J[53] + x86 * J[54];
    double j53 = x79 * J[52] + x83 * J[53] + x87 * J[54];
    double j54 = x81 * J[52] + x85 * J[53] + x88 * J[54];
    double j55 = x89 * J[52] + x90 * J[53] + x91 * J[54];

    J[4] = j4;
    J[5] = j5;
    J[6] = j6;
    J[7] += j7;

    J[12] = j12;
    J[13] = j13;
    J[14] = j14;
    J[15] += j15;

    J[20] = j20;
    J[21] = j21;
    J[22] = j22;
    J[23] += j23;

    J[31] += h * lambda * std::pow(m, 2) / dtds;

    J[36] = j36;
    J[37] = j37;
    J[38] = j38;
    J[39] += j39;

    J[44] = j44;
    J[45] = j45;
    J[46] = j46;
    J[47] += j47;

    J[52] = j52;
    J[53] = j53;
    J[54] = j54;
    J[55] += j55;
  }

  return true;
}

}  // namespace

SympyStepper::SympyStepper(std::shared_ptr<const MagneticFieldProvider> bField,
                           double overstepLimit)
    : m_bField(std::move(bField)), m_overstepLimit(overstepLimit) {}

SympyStepper::State SympyStepper::makeState(
    std::reference_wrapper<const GeometryContext> gctx,
    std::reference_wrapper<const MagneticFieldContext> mctx,
    const BoundTrackParameters& par, double ssize) const {
  return State{gctx, m_bField->makeCache(mctx), par, ssize};
}

void SympyStepper::resetState(State& state, const BoundVector& boundParams,
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
SympyStepper::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  return detail::boundState(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated,
      freeToBoundCorrection);
}

std::tuple<CurvilinearTrackParameters, BoundMatrix, double>
SympyStepper::curvilinearState(State& state, bool transportCov) const {
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated);
}

void SympyStepper::update(State& state, const FreeVector& freeParams,
                          const BoundVector& /*boundParams*/,
                          const Covariance& covariance,
                          const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
}

void SympyStepper::update(State& state, const Vector3& uposition,
                          const Vector3& udirection, double qop,
                          double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qop;
}

void SympyStepper::transportCovarianceToCurvilinear(State& state) const {
  detail::transportCovarianceToCurvilinear(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars.template segment<3>(eFreeDir0));
}

void SympyStepper::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::transportCovarianceToBound(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, freeToBoundCorrection);
}

Result<double> SympyStepper::stepImpl(
    State& state, Direction stepDirection, double stepTolerance,
    double stepSizeCutOff, std::size_t maxRungeKuttaStepTrials) const {
  // Runge-Kutta integrator state
  auto& sd = state.stepData;

  auto pos = position(state);
  auto dir = direction(state);
  double qop = qOverP(state);
  double m = particleHypothesis(state).mass();
  double p = absoluteMomentum(state);

  auto getB = [&](const double* p) -> Vector3 {
    auto fieldRes = getField(state, {p[0], p[1], p[2]});
    return *fieldRes;
  };

  double h = state.stepSize.value() * stepDirection;
  std::size_t nStepTrials = 0;
  double error_estimate = 0.;

  while (true) {
    bool ok = rk4(pos.data(), dir.data(), h, qop, m, p,
                  state.pars.template segment<3>(eFreePos0).data(),
                  state.pars.template segment<3>(eFreeDir0).data(),
                  state.pars.template segment<1>(eFreeTime).data(),
                  state.jacTransport.data(), getB, state.covTransport);

    if (ok) {
      break;
    }

    /*
    // double std::sqrt is 3x faster than std::pow
    const double stepSizeScaling = std::clamp(
        std::sqrt(std::sqrt(stepTolerance / std::abs(2. * error_estimate))),
        0.25, 4.0);
    h *= stepSizeScaling;
    */

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

  return h;
}

void SympyStepper::setIdentityJacobian(State& state) const {
  state.jacobian = BoundMatrix::Identity();
}

}  // namespace Acts
