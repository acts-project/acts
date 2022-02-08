// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"

template <typename E, typename A>
Acts::EigenStepper<E, A>::EigenStepper(
    std::shared_ptr<const MagneticFieldProvider> bField)
    : m_bField(std::move(bField)) {}

template <typename E, typename A>
template <typename charge_t>
auto Acts::EigenStepper<E, A>::makeState(
    std::reference_wrapper<const GeometryContext> gctx,
    std::reference_wrapper<const MagneticFieldContext> mctx,
    const SingleBoundTrackParameters<charge_t>& par, NavigationDirection ndir,
    double ssize, double stolerance) const -> State {
  return State{gctx, m_bField->makeCache(mctx), par, ndir, ssize, stolerance};
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::resetState(State& state,
                                          const BoundVector& boundParams,
                                          const BoundSymMatrix& cov,
                                          const Surface& surface,
                                          const NavigationDirection navDir,
                                          const double stepSize) const {
  // Update the stepping state
  update(state,
         detail::transformBoundToFreeParameters(surface, state.geoContext,
                                                boundParams),
         boundParams, cov, surface);
  state.navDir = navDir;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  state.jacToGlobal =
      surface.boundToFreeJacobian(state.geoContext, boundParams);
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

template <typename E, typename A>
auto Acts::EigenStepper<E, A>::boundState(State& state, const Surface& surface,
                                          bool transportCov) const
    -> Result<BoundState> {
  return detail::boundState(
      state.geoContext, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars,
      state.covTransport && transportCov, state.pathAccumulated, surface);
}

template <typename E, typename A>
auto Acts::EigenStepper<E, A>::curvilinearState(State& state,
                                                bool transportCov) const
    -> CurvilinearState {
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars, state.covTransport && transportCov,
      state.pathAccumulated);
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::update(State& state,
                                      const FreeVector& freeParams,
                                      const BoundVector& boundParams,
                                      const Covariance& covariance,
                                      const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal =
      surface.boundToFreeJacobian(state.geoContext, boundParams);
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::update(State& state, const Vector3& uposition,
                                      const Vector3& udirection, double up,
                                      double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = (state.q != 0. ? state.q / up : 1. / up);
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::transportCovarianceToCurvilinear(
    State& state) const {
  detail::transportCovarianceToCurvilinear(state.cov, state.jacobian,
                                           state.jacTransport, state.derivative,
                                           state.jacToGlobal, direction(state));
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::transportCovarianceToBound(
    State& state, const Surface& surface) const {
  detail::transportCovarianceToBound(
      state.geoContext.get(), state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, surface);
}

inline void mult_8x8_8x8(double* D_result, const double* A, const double* B) {
  D_result[0] = A[0] * B[0] + A[1] * B[8] + A[2] * B[16] + A[3] * B[24] +
                A[4] * B[32] + A[5] * B[40] + A[6] * B[48] + A[7] * B[56];
  D_result[1] = A[0] * B[1] + A[1] * B[9] + A[2] * B[17] + A[3] * B[25] +
                A[4] * B[33] + A[5] * B[41] + A[6] * B[49] + A[7] * B[57];
  D_result[2] = A[0] * B[2] + A[1] * B[10] + A[2] * B[18] + A[3] * B[26] +
                A[4] * B[34] + A[5] * B[42] + A[6] * B[50] + A[7] * B[58];
  D_result[3] = A[0] * B[3] + A[1] * B[11] + A[2] * B[19] + A[3] * B[27] +
                A[4] * B[35] + A[5] * B[43] + A[6] * B[51] + A[7] * B[59];
  D_result[4] = A[0] * B[4] + A[1] * B[12] + A[2] * B[20] + A[3] * B[28] +
                A[4] * B[36] + A[5] * B[44] + A[6] * B[52] + A[7] * B[60];
  D_result[5] = A[0] * B[5] + A[1] * B[13] + A[2] * B[21] + A[3] * B[29] +
                A[4] * B[37] + A[5] * B[45] + A[6] * B[53] + A[7] * B[61];
  D_result[6] = A[0] * B[6] + A[1] * B[14] + A[2] * B[22] + A[3] * B[30] +
                A[4] * B[38] + A[5] * B[46] + A[6] * B[54] + A[7] * B[62];
  D_result[7] = A[0] * B[7] + A[1] * B[15] + A[2] * B[23] + A[3] * B[31] +
                A[4] * B[39] + A[5] * B[47] + A[6] * B[55] + A[7] * B[63];
  D_result[8] = A[8] * B[0] + A[9] * B[8] + A[10] * B[16] + A[11] * B[24] +
                A[12] * B[32] + A[13] * B[40] + A[14] * B[48] + A[15] * B[56];
  D_result[9] = A[8] * B[1] + A[9] * B[9] + A[10] * B[17] + A[11] * B[25] +
                A[12] * B[33] + A[13] * B[41] + A[14] * B[49] + A[15] * B[57];
  D_result[10] = A[8] * B[2] + A[9] * B[10] + A[10] * B[18] + A[11] * B[26] +
                 A[12] * B[34] + A[13] * B[42] + A[14] * B[50] + A[15] * B[58];
  D_result[11] = A[8] * B[3] + A[9] * B[11] + A[10] * B[19] + A[11] * B[27] +
                 A[12] * B[35] + A[13] * B[43] + A[14] * B[51] + A[15] * B[59];
  D_result[12] = A[8] * B[4] + A[9] * B[12] + A[10] * B[20] + A[11] * B[28] +
                 A[12] * B[36] + A[13] * B[44] + A[14] * B[52] + A[15] * B[60];
  D_result[13] = A[8] * B[5] + A[9] * B[13] + A[10] * B[21] + A[11] * B[29] +
                 A[12] * B[37] + A[13] * B[45] + A[14] * B[53] + A[15] * B[61];
  D_result[14] = A[8] * B[6] + A[9] * B[14] + A[10] * B[22] + A[11] * B[30] +
                 A[12] * B[38] + A[13] * B[46] + A[14] * B[54] + A[15] * B[62];
  D_result[15] = A[8] * B[7] + A[9] * B[15] + A[10] * B[23] + A[11] * B[31] +
                 A[12] * B[39] + A[13] * B[47] + A[14] * B[55] + A[15] * B[63];
  D_result[16] = A[16] * B[0] + A[17] * B[8] + A[18] * B[16] + A[19] * B[24] +
                 A[20] * B[32] + A[21] * B[40] + A[22] * B[48] + A[23] * B[56];
  D_result[17] = A[16] * B[1] + A[17] * B[9] + A[18] * B[17] + A[19] * B[25] +
                 A[20] * B[33] + A[21] * B[41] + A[22] * B[49] + A[23] * B[57];
  D_result[18] = A[16] * B[2] + A[17] * B[10] + A[18] * B[18] + A[19] * B[26] +
                 A[20] * B[34] + A[21] * B[42] + A[22] * B[50] + A[23] * B[58];
  D_result[19] = A[16] * B[3] + A[17] * B[11] + A[18] * B[19] + A[19] * B[27] +
                 A[20] * B[35] + A[21] * B[43] + A[22] * B[51] + A[23] * B[59];
  D_result[20] = A[16] * B[4] + A[17] * B[12] + A[18] * B[20] + A[19] * B[28] +
                 A[20] * B[36] + A[21] * B[44] + A[22] * B[52] + A[23] * B[60];
  D_result[21] = A[16] * B[5] + A[17] * B[13] + A[18] * B[21] + A[19] * B[29] +
                 A[20] * B[37] + A[21] * B[45] + A[22] * B[53] + A[23] * B[61];
  D_result[22] = A[16] * B[6] + A[17] * B[14] + A[18] * B[22] + A[19] * B[30] +
                 A[20] * B[38] + A[21] * B[46] + A[22] * B[54] + A[23] * B[62];
  D_result[23] = A[16] * B[7] + A[17] * B[15] + A[18] * B[23] + A[19] * B[31] +
                 A[20] * B[39] + A[21] * B[47] + A[22] * B[55] + A[23] * B[63];
  D_result[24] = A[24] * B[0] + A[25] * B[8] + A[26] * B[16] + A[27] * B[24] +
                 A[28] * B[32] + A[29] * B[40] + A[30] * B[48] + A[31] * B[56];
  D_result[25] = A[24] * B[1] + A[25] * B[9] + A[26] * B[17] + A[27] * B[25] +
                 A[28] * B[33] + A[29] * B[41] + A[30] * B[49] + A[31] * B[57];
  D_result[26] = A[24] * B[2] + A[25] * B[10] + A[26] * B[18] + A[27] * B[26] +
                 A[28] * B[34] + A[29] * B[42] + A[30] * B[50] + A[31] * B[58];
  D_result[27] = A[24] * B[3] + A[25] * B[11] + A[26] * B[19] + A[27] * B[27] +
                 A[28] * B[35] + A[29] * B[43] + A[30] * B[51] + A[31] * B[59];
  D_result[28] = A[24] * B[4] + A[25] * B[12] + A[26] * B[20] + A[27] * B[28] +
                 A[28] * B[36] + A[29] * B[44] + A[30] * B[52] + A[31] * B[60];
  D_result[29] = A[24] * B[5] + A[25] * B[13] + A[26] * B[21] + A[27] * B[29] +
                 A[28] * B[37] + A[29] * B[45] + A[30] * B[53] + A[31] * B[61];
  D_result[30] = A[24] * B[6] + A[25] * B[14] + A[26] * B[22] + A[27] * B[30] +
                 A[28] * B[38] + A[29] * B[46] + A[30] * B[54] + A[31] * B[62];
  D_result[31] = A[24] * B[7] + A[25] * B[15] + A[26] * B[23] + A[27] * B[31] +
                 A[28] * B[39] + A[29] * B[47] + A[30] * B[55] + A[31] * B[63];
  D_result[32] = A[32] * B[0] + A[33] * B[8] + A[34] * B[16] + A[35] * B[24] +
                 A[36] * B[32] + A[37] * B[40] + A[38] * B[48] + A[39] * B[56];
  D_result[33] = A[32] * B[1] + A[33] * B[9] + A[34] * B[17] + A[35] * B[25] +
                 A[36] * B[33] + A[37] * B[41] + A[38] * B[49] + A[39] * B[57];
  D_result[34] = A[32] * B[2] + A[33] * B[10] + A[34] * B[18] + A[35] * B[26] +
                 A[36] * B[34] + A[37] * B[42] + A[38] * B[50] + A[39] * B[58];
  D_result[35] = A[32] * B[3] + A[33] * B[11] + A[34] * B[19] + A[35] * B[27] +
                 A[36] * B[35] + A[37] * B[43] + A[38] * B[51] + A[39] * B[59];
  D_result[36] = A[32] * B[4] + A[33] * B[12] + A[34] * B[20] + A[35] * B[28] +
                 A[36] * B[36] + A[37] * B[44] + A[38] * B[52] + A[39] * B[60];
  D_result[37] = A[32] * B[5] + A[33] * B[13] + A[34] * B[21] + A[35] * B[29] +
                 A[36] * B[37] + A[37] * B[45] + A[38] * B[53] + A[39] * B[61];
  D_result[38] = A[32] * B[6] + A[33] * B[14] + A[34] * B[22] + A[35] * B[30] +
                 A[36] * B[38] + A[37] * B[46] + A[38] * B[54] + A[39] * B[62];
  D_result[39] = A[32] * B[7] + A[33] * B[15] + A[34] * B[23] + A[35] * B[31] +
                 A[36] * B[39] + A[37] * B[47] + A[38] * B[55] + A[39] * B[63];
  D_result[40] = A[40] * B[0] + A[41] * B[8] + A[42] * B[16] + A[43] * B[24] +
                 A[44] * B[32] + A[45] * B[40] + A[46] * B[48] + A[47] * B[56];
  D_result[41] = A[40] * B[1] + A[41] * B[9] + A[42] * B[17] + A[43] * B[25] +
                 A[44] * B[33] + A[45] * B[41] + A[46] * B[49] + A[47] * B[57];
  D_result[42] = A[40] * B[2] + A[41] * B[10] + A[42] * B[18] + A[43] * B[26] +
                 A[44] * B[34] + A[45] * B[42] + A[46] * B[50] + A[47] * B[58];
  D_result[43] = A[40] * B[3] + A[41] * B[11] + A[42] * B[19] + A[43] * B[27] +
                 A[44] * B[35] + A[45] * B[43] + A[46] * B[51] + A[47] * B[59];
  D_result[44] = A[40] * B[4] + A[41] * B[12] + A[42] * B[20] + A[43] * B[28] +
                 A[44] * B[36] + A[45] * B[44] + A[46] * B[52] + A[47] * B[60];
  D_result[45] = A[40] * B[5] + A[41] * B[13] + A[42] * B[21] + A[43] * B[29] +
                 A[44] * B[37] + A[45] * B[45] + A[46] * B[53] + A[47] * B[61];
  D_result[46] = A[40] * B[6] + A[41] * B[14] + A[42] * B[22] + A[43] * B[30] +
                 A[44] * B[38] + A[45] * B[46] + A[46] * B[54] + A[47] * B[62];
  D_result[47] = A[40] * B[7] + A[41] * B[15] + A[42] * B[23] + A[43] * B[31] +
                 A[44] * B[39] + A[45] * B[47] + A[46] * B[55] + A[47] * B[63];
  D_result[48] = A[48] * B[0] + A[49] * B[8] + A[50] * B[16] + A[51] * B[24] +
                 A[52] * B[32] + A[53] * B[40] + A[54] * B[48] + A[55] * B[56];
  D_result[49] = A[48] * B[1] + A[49] * B[9] + A[50] * B[17] + A[51] * B[25] +
                 A[52] * B[33] + A[53] * B[41] + A[54] * B[49] + A[55] * B[57];
  D_result[50] = A[48] * B[2] + A[49] * B[10] + A[50] * B[18] + A[51] * B[26] +
                 A[52] * B[34] + A[53] * B[42] + A[54] * B[50] + A[55] * B[58];
  D_result[51] = A[48] * B[3] + A[49] * B[11] + A[50] * B[19] + A[51] * B[27] +
                 A[52] * B[35] + A[53] * B[43] + A[54] * B[51] + A[55] * B[59];
  D_result[52] = A[48] * B[4] + A[49] * B[12] + A[50] * B[20] + A[51] * B[28] +
                 A[52] * B[36] + A[53] * B[44] + A[54] * B[52] + A[55] * B[60];
  D_result[53] = A[48] * B[5] + A[49] * B[13] + A[50] * B[21] + A[51] * B[29] +
                 A[52] * B[37] + A[53] * B[45] + A[54] * B[53] + A[55] * B[61];
  D_result[54] = A[48] * B[6] + A[49] * B[14] + A[50] * B[22] + A[51] * B[30] +
                 A[52] * B[38] + A[53] * B[46] + A[54] * B[54] + A[55] * B[62];
  D_result[55] = A[48] * B[7] + A[49] * B[15] + A[50] * B[23] + A[51] * B[31] +
                 A[52] * B[39] + A[53] * B[47] + A[54] * B[55] + A[55] * B[63];
  D_result[56] = A[56] * B[0] + A[57] * B[8] + A[58] * B[16] + A[59] * B[24] +
                 A[60] * B[32] + A[61] * B[40] + A[62] * B[48] + A[63] * B[56];
  D_result[57] = A[56] * B[1] + A[57] * B[9] + A[58] * B[17] + A[59] * B[25] +
                 A[60] * B[33] + A[61] * B[41] + A[62] * B[49] + A[63] * B[57];
  D_result[58] = A[56] * B[2] + A[57] * B[10] + A[58] * B[18] + A[59] * B[26] +
                 A[60] * B[34] + A[61] * B[42] + A[62] * B[50] + A[63] * B[58];
  D_result[59] = A[56] * B[3] + A[57] * B[11] + A[58] * B[19] + A[59] * B[27] +
                 A[60] * B[35] + A[61] * B[43] + A[62] * B[51] + A[63] * B[59];
  D_result[60] = A[56] * B[4] + A[57] * B[12] + A[58] * B[20] + A[59] * B[28] +
                 A[60] * B[36] + A[61] * B[44] + A[62] * B[52] + A[63] * B[60];
  D_result[61] = A[56] * B[5] + A[57] * B[13] + A[58] * B[21] + A[59] * B[29] +
                 A[60] * B[37] + A[61] * B[45] + A[62] * B[53] + A[63] * B[61];
  D_result[62] = A[56] * B[6] + A[57] * B[14] + A[58] * B[22] + A[59] * B[30] +
                 A[60] * B[38] + A[61] * B[46] + A[62] * B[54] + A[63] * B[62];
  D_result[63] = A[56] * B[7] + A[57] * B[15] + A[58] * B[23] + A[59] * B[31] +
                 A[60] * B[39] + A[61] * B[47] + A[62] * B[55] + A[63] * B[63];
}

template <typename E, typename A>
template <typename propagator_state_t>
Acts::Result<double> Acts::EigenStepper<E, A>::step(
    propagator_state_t& state) const {
  using namespace UnitLiterals;

  // Runge-Kutta integrator state
  auto& sd = state.stepping.stepData;
  double error_estimate = 0.;
  double h2, half_h;

  auto pos = position(state.stepping);
  auto dir = direction(state.stepping);

  // First Runge-Kutta point (at current position)
  auto fieldRes = getField(state.stepping, pos);
  if (!fieldRes.ok()) {
    return fieldRes.error();
  }
  sd.B_first = *fieldRes;
  if (!state.stepping.extension.validExtensionForStep(state, *this) ||
      !state.stepping.extension.k1(state, *this, sd.k1, sd.B_first, sd.kQoP)) {
    return 0.;
  }

  // The following functor starts to perform a Runge-Kutta step of a certain
  // size, going up to the point where it can return an estimate of the local
  // integration error. The results are stated in the local variables above,
  // allowing integration to continue once the error is deemed satisfactory
  const auto tryRungeKuttaStep = [&](const ConstrainedStep& h) -> Result<bool> {
    // helpers because bool and std::error_code are ambiguous
    constexpr auto success = &Result<bool>::success;
    constexpr auto failure = &Result<bool>::failure;

    // State the square and half of the step size
    h2 = h * h;
    half_h = h * 0.5;

    // Second Runge-Kutta point
    const Vector3 pos1 = pos + half_h * dir + h2 * 0.125 * sd.k1;
    auto field = getField(state.stepping, pos1);
    if (!field.ok()) {
      return failure(field.error());
    }
    sd.B_middle = *field;

    if (!state.stepping.extension.k2(state, *this, sd.k2, sd.B_middle, sd.kQoP,
                                     half_h, sd.k1)) {
      return success(false);
    }

    // Third Runge-Kutta point
    if (!state.stepping.extension.k3(state, *this, sd.k3, sd.B_middle, sd.kQoP,
                                     half_h, sd.k2)) {
      return success(false);
    }

    // Last Runge-Kutta point
    const Vector3 pos2 = pos + h * dir + h2 * 0.5 * sd.k3;
    field = getField(state.stepping, pos2);
    if (!field.ok()) {
      return failure(field.error());
    }
    sd.B_last = *field;
    if (!state.stepping.extension.k4(state, *this, sd.k4, sd.B_last, sd.kQoP, h,
                                     sd.k3)) {
      return success(false);
    }

    // Compute and check the local integration error estimate
    error_estimate = std::max(
        h2 * ((sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>() +
              std::abs(sd.kQoP[0] - sd.kQoP[1] - sd.kQoP[2] + sd.kQoP[3])),
        1e-20);

    return success(error_estimate <= state.options.tolerance);
  };

  double stepSizeScaling = 1.;
  size_t nStepTrials = 0;
  // Select and adjust the appropriate Runge-Kutta step size as given
  // ATL-SOFT-PUB-2009-001
  while (true) {
    auto res = tryRungeKuttaStep(state.stepping.stepSize);
    if (!res.ok()) {
      return res.error();
    }
    if (!!res.value()) {
      break;
    }

    stepSizeScaling = std::min(
        std::max(0.25, std::sqrt(std::sqrt((state.options.tolerance /
                                            std::abs(2. * error_estimate))))),
        4.);

    state.stepping.stepSize = state.stepping.stepSize * stepSizeScaling;

    // If step size becomes too small the particle remains at the initial
    // place
    if (std::abs(state.stepping.stepSize) <
        std::abs(state.options.stepSizeCutOff)) {
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

    // // for moment, only update the transport part
    // // state.stepping.jacTransport = D * state.stepping.jacTransport;
    // //
    Acts::FreeMatrix jac = state.stepping.jacTransport;
    mult_8x8_8x8(state.stepping.jacTransport.data(), D.data(), jac.data());
  } else {
    if (!state.stepping.extension.finalize(state, *this, h)) {
      return EigenStepperError::StepInvalid;
    }
  }

  // Update the track parameters according to the equations of motion
  state.stepping.pars.template segment<3>(eFreePos0) +=
      h * dir + h2 / 6. * (sd.k1 + sd.k2 + sd.k3);
  state.stepping.pars.template segment<3>(eFreeDir0) +=
      h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
  (state.stepping.pars.template segment<3>(eFreeDir0)).normalize();

  if (state.stepping.covTransport) {
    state.stepping.derivative.template head<3>() =
        state.stepping.pars.template segment<3>(eFreeDir0);
    state.stepping.derivative.template segment<3>(4) = sd.k4;
  }
  state.stepping.pathAccumulated += h;
  if (state.stepping.stepSize.currentType() ==
      ConstrainedStep::Type::accuracy) {
    state.stepping.stepSize =
        state.stepping.stepSize *
        std::min(
            std::max(0.25, std::sqrt(std::sqrt((state.options.tolerance /
                                                std::abs(error_estimate))))),
            4.);
  }
  return h;
}
