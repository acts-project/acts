// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/HelixStepper.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>

namespace Acts {

namespace {

/// Even functions of the turning angle theta = |omega| * |h| that the helix
/// and its jacobian are built from. Each one is normalised so that it is
/// finite at theta = 0.
struct HelixCoefficients {
  /// sin(theta) / theta
  double f1 = 0;
  /// (1 - cos(theta)) / theta^2
  double f2 = 0;
  /// (theta - sin(theta)) / theta^3
  double g2 = 0;
  /// (sin(theta) - theta * cos(theta)) / theta^3
  double h1 = 0;
  /// (theta^2 / 2 + 1 - cos(theta) - theta * sin(theta)) / theta^4
  double h2 = 0;
};

HelixCoefficients helixCoefficients(double theta2) {
  HelixCoefficients c;
  // The closed forms cancel for a small angle, so use the series there. The
  // first omitted term is below 1e-17 at the threshold.
  if (theta2 < 1e-3) {
    const double t2 = theta2;
    const double t4 = t2 * t2;
    const double t6 = t4 * t2;
    c.f1 = 1. - t2 / 6. + t4 / 120. - t6 / 5040.;
    c.f2 = 1. / 2. - t2 / 24. + t4 / 720. - t6 / 40320.;
    c.g2 = 1. / 6. - t2 / 120. + t4 / 5040. - t6 / 362880.;
    c.h1 = 1. / 3. - t2 / 30. + t4 / 840. - t6 / 45360.;
    c.h2 = 1. / 8. - t2 / 144. + t4 / 5760. - t6 / 403200.;
    return c;
  }
  const double theta = std::sqrt(theta2);
  const double sinTheta = std::sin(theta);
  const double cosTheta = std::cos(theta);
  const double halfSin = std::sin(0.5 * theta);
  // 1 - cos(theta) without the cancellation
  const double oneMinusCos = 2. * halfSin * halfSin;
  const double theta3 = theta2 * theta;
  c.f1 = sinTheta / theta;
  c.f2 = oneMinusCos / theta2;
  c.g2 = (theta - sinTheta) / theta3;
  c.h1 = (sinTheta - theta * cosTheta) / theta3;
  c.h2 = (0.5 * theta2 + oneMinusCos - theta * sinTheta) / (theta2 * theta2);
  return c;
}

/// Cross product matrix W with W * v = omega x v
SquareMatrix3 crossMatrix(const Vector3& omega) {
  SquareMatrix3 w;
  // clang-format off
  w <<         0., -omega.z(),  omega.y(),
        omega.z(),         0., -omega.x(),
       -omega.y(),  omega.x(),         0.;
  // clang-format on
  return w;
}

}  // namespace

HelixStepper::HelixStepper(std::shared_ptr<const MagneticFieldProvider> bField)
    : m_bField(std::move(bField)) {}

HelixStepper::HelixStepper(const Config& config) : m_bField(config.bField) {}

HelixStepper::State HelixStepper::makeState(const Options& options) const {
  State state{options, m_bField->makeCache(options.magFieldContext)};
  return state;
}

void HelixStepper::initialize(State& state, const BoundParameters& par) const {
  initialize(state, par.parameters(), par.covariance(),
             par.particleHypothesis(), par.referenceSurface());
}

void HelixStepper::initialize(State& state, const BoundVector& boundParams,
                              const std::optional<BoundMatrix>& cov,
                              ParticleHypothesis particleHypothesis,
                              const Surface& surface) const {
  FreeVector freeParams = transformBoundToFreeParameters(
      surface, state.options.geoContext, boundParams);

  state.particleHypothesis = particleHypothesis;

  state.pathAccumulated = 0;
  state.nSteps = 0;
  state.stepSize = ConstrainedStep();
  state.stepSize.setAccuracy(state.options.initialStepSize);
  state.stepSize.setUser(state.options.maxStepSize);
  state.previousStepSize = 0;
  state.statistics = StepperStatistics();

  state.pars = freeParams;
  state.field.reset();

  state.cov = cov;
  if (state.cov.has_value()) {
    state.jacToGlobal = surface.boundToFreeJacobian(
        state.options.geoContext, freeParams.segment<3>(eFreePos0),
        freeParams.segment<3>(eFreeDir0));
    state.jacTransport = FreeMatrix::Identity();
    state.derivative = FreeVector::Zero();
  }
}

Result<HelixStepper::BoundParameters> HelixStepper::boundParameters(
    const State& state, const Surface& surface) const {
  return detail::boundParameters(state.options.geoContext, surface, state.pars,
                                 state.cov, state.particleHypothesis);
}

bool HelixStepper::prepareCurvilinearState(State& state) const {
  // The derivatives are only missing if no step was executed yet.
  if (state.pathAccumulated != 0) {
    return true;
  }

  if (!state.field.has_value()) {
    auto fieldRes = getField(state, position(state));
    if (!fieldRes.ok()) {
      return false;
    }
    state.field = *fieldRes;
  }

  const double qop = qOverP(state);
  const double mass = state.particleHypothesis.mass();
  state.derivative.segment<3>(eFreePos0) = direction(state);
  state.derivative[eFreeTime] = fastHypot(1., mass / absoluteMomentum(state));
  state.derivative.segment<3>(eFreeDir0) =
      qop * direction(state).cross(*state.field);
  state.derivative[eFreeQOverP] = 0.;
  return true;
}

HelixStepper::BoundParameters HelixStepper::curvilinearParameters(
    const State& state) const {
  return detail::curvilinearParameters(state.pars, state.cov,
                                       state.particleHypothesis);
}

void HelixStepper::update(State& state, const FreeVector& freeParams,
                          const BoundVector& /*boundParams*/,
                          const Covariance& covariance,
                          const Surface& surface) const {
  state.pars = freeParams;
  state.field.reset();
  if (state.cov.has_value()) {
    state.cov = covariance;
  }
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.options.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

void HelixStepper::update(State& state, const Vector3& uposition,
                          const Vector3& udirection, double qop,
                          double time) const {
  // Material interactions keep the position, so the cached field stays valid.
  if (uposition != position(state)) {
    state.field.reset();
  }
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qop;
}

HelixStepper::Jacobian HelixStepper::transportToCurvilinear(
    State& state) const {
  Jacobian jacobian = Jacobian::Identity();
  if (!state.cov.has_value()) {
    return jacobian;
  }
  detail::transportCovarianceToCurvilinear(
      *state.cov, jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, std::nullopt, direction(state));
  return jacobian;
}

Result<HelixStepper::Jacobian> HelixStepper::transportToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  if (!surface.isOnSurface(state.options.geoContext, position(state),
                           direction(state), BoundaryTolerance::Infinite())) {
    return Result<Jacobian>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }

  Jacobian jacobian = Jacobian::Identity();
  if (!state.cov.has_value()) {
    return Result<Jacobian>::success(jacobian);
  }
  detail::transportCovarianceToBound(
      state.options.geoContext, surface, *state.cov, jacobian,
      state.jacTransport, state.derivative, state.jacToGlobal, std::nullopt,
      state.pars, freeToBoundCorrection);
  return Result<Jacobian>::success(jacobian);
}

Result<double> HelixStepper::step(State& state, Direction propDir,
                                  const IVolumeMaterial* material) const {
  static_cast<void>(material);

  const Vector3 pos = position(state);
  const Vector3 dir = direction(state);
  const double qop = qOverP(state);

  if (!state.field.has_value()) {
    auto fieldRes = getField(state, pos);
    if (!fieldRes.ok()) {
      return fieldRes.error();
    }
    state.field = *fieldRes;
  }
  const Vector3 startField = *state.field;

  // The equation of motion is dT/ds = T x omega. Its solution is a rotation
  // of T about omega:
  //   T(h) = T - h f1 a + h^2 f2 b
  //   r(h) = r + h T - h^2 f2 a + h^3 g2 b
  // with a = omega x T and b = omega x a.
  const Vector3 omega = qop * startField;
  const Vector3 a = omega.cross(dir);
  const Vector3 b = omega.cross(a);
  const double omega2 = omega.squaredNorm();

  const auto calcStepSizeScaling = [&](const double errorEstimate) -> double {
    constexpr double lower = 0.25;
    constexpr double upper = 4.0;
    // The error of a linear field change along the step grows with h^3
    const double x = std::cbrt(state.options.stepTolerance / errorEstimate);
    return std::clamp(x, lower, upper);
  };

  const double initialH = state.stepSize.value() * propDir;
  double h = initialH;
  HelixCoefficients coeffs;
  Vector3 endPos;
  Vector3 endDir;
  Vector3 endField;
  double errorEstimate = 0.;
  std::size_t nStepTrials = 0;

  while (true) {
    ++nStepTrials;
    ++state.statistics.nAttemptedSteps;

    coeffs = helixCoefficients(omega2 * h * h);
    const double h2 = h * h;
    endDir = dir - h * coeffs.f1 * a + h2 * coeffs.f2 * b;
    endPos = pos + h * dir - h2 * coeffs.f2 * a + h2 * h * coeffs.g2 * b;

    // The end field is the start field of the next step, so an accepted step
    // costs no extra lookup.
    auto fieldRes = getField(state, endPos);
    if (!fieldRes.ok()) {
      return fieldRes.error();
    }
    endField = *fieldRes;

    // Position deviation from a field that changes linearly along the step
    errorEstimate = std::max(
        1e-20, h2 / 6. * (qop * dir.cross(endField - startField)).norm());

    if (errorEstimate <= state.options.stepTolerance) {
      break;
    }

    ++state.statistics.nRejectedSteps;

    h *= calcStepSizeScaling(errorEstimate);

    if (std::abs(h) < std::abs(state.options.stepSizeCutOff)) {
      return EigenStepperError::StepSizeStalled;
    }

    if (nStepTrials > state.options.maxRungeKuttaStepTrials) {
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
  }

  const double mass = state.particleHypothesis.mass();
  const double absMom = absoluteMomentum(state);
  // time propagates along distance as 1/b = sqrt(1 + m²/p²)
  const double dtds = fastHypot(1., mass / absMom);

  if (state.cov.has_value()) {
    const double h2 = h * h;
    const SquareMatrix3 w = crossMatrix(omega);
    const SquareMatrix3 w2 = w * w;

    // The helix does not depend on the start position and time, so the step
    // jacobian D has D11 = I and D21 = 0. Only the right half is non-trivial.
    Matrix<4, 4> dTopRight = Matrix<4, 4>::Zero();
    Matrix<4, 4> dBottomRight = Matrix<4, 4>::Identity();

    // dr/dT and dT/dT
    dTopRight.topLeftCorner<3, 3>() = h * SquareMatrix3::Identity() -
                                      h2 * coeffs.f2 * w +
                                      h2 * h * coeffs.g2 * w2;
    dBottomRight.topLeftCorner<3, 3>() =
        SquareMatrix3::Identity() - h * coeffs.f1 * w + h2 * coeffs.f2 * w2;

    // omega = qop * B, so d/d(qop) = B x d/d(omega) applied to the
    // integrals of T:
    //   dT(h)/d(qop) = h T(h) x B
    //   dr(h)/d(qop) = (h^2/2 T - h^3 h1 a + h^4 h2 b) x B
    dTopRight.block<3, 1>(0, 3) =
        (0.5 * h2 * dir - h2 * h * coeffs.h1 * a + h2 * h2 * coeffs.h2 * b)
            .cross(startField);
    dBottomRight.block<3, 1>(0, 3) = h * endDir.cross(startField);

    // dt/d(qop) = h m^2 (q/p) / (q^2 dt/ds). A neutral particle has a fixed
    // momentum hypothesis, so its time does not depend on qop.
    const double q = charge(state);
    dTopRight(3, 3) = q != 0. ? h * mass * mass * qop / (q * q * dtds) : 0.;

    // See EigenStepper::step for the blocked multiplication.
    state.jacTransport.template topRightCorner<4, 4>() +=
        dTopRight * state.jacTransport.template bottomRightCorner<4, 4>();
    state.jacTransport.template bottomRightCorner<4, 4>() =
        (dBottomRight * state.jacTransport.template bottomRightCorner<4, 4>())
            .eval();

    state.derivative.template segment<3>(eFreePos0) = endDir;
    state.derivative[eFreeTime] = dtds;
    state.derivative.template segment<3>(eFreeDir0) = endDir.cross(omega);
    state.derivative[eFreeQOverP] = 0.;
  }

  state.pars.template segment<3>(eFreePos0) = endPos;
  state.pars.template segment<3>(eFreeDir0) = endDir;
  state.pars[eFreeTime] += h * dtds;
  state.field = endField;

  state.pathAccumulated += h;
  ++state.nSteps;

  ++state.statistics.nSuccessfulSteps;
  if (propDir != Direction::fromScalarZeroAsPositive(initialH)) {
    ++state.statistics.nReverseSteps;
  }
  state.statistics.pathLength += h;
  state.statistics.absolutePathLength += std::abs(h);

  const double nextAccuracy = std::abs(h * calcStepSizeScaling(errorEstimate));
  const double previousAccuracy = std::abs(state.stepSize.accuracy());
  const double initialStepLength = std::abs(initialH);
  if (nextAccuracy < initialStepLength || nextAccuracy > previousAccuracy) {
    state.stepSize.setAccuracy(nextAccuracy);
  }

  return h;
}

}  // namespace Acts
