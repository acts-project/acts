// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <numbers>
#include <optional>
#include <utility>

namespace Acts::Experimental {

/// Space point as a space point filter reads it.
struct SpacePointFilterInput final {
  float x{};
  float y{};
  float z{};
  /// Transverse distance from the beamline.
  float r{};
  /// Whether r is measured at a known z (disc), instead of z at a known r
  /// (cylinder).
  bool measuresR{false};
  /// Variance across the track in the transverse plane.
  float varianceTransverse{};
  /// Variance of the measured longitudinal coordinate, z or r.
  float varianceLongitudinal{};
  /// Half of the outer strip for a strip space point; the crossing is
  /// unconstrained along it.
  std::optional<std::array<float, 3>> stripHalfVector;
  /// Radial field at the space point. The field is taken axially symmetric,
  /// without an azimuthal component, as in a solenoid.
  float br{};
  /// Field along z at the space point.
  float bz{};
  /// Material at the space point in radiation lengths, at normal incidence.
  float thicknessInX0{};
};

/// Track parameters at the last space point, along the fit direction. In
/// float like the filters; the caller widens.
struct SpacePointFilterParameters final {
  /// Position at the last space point.
  Eigen::Vector3f position{Eigen::Vector3f::Zero()};
  /// Azimuthal angle of the direction.
  float phi{0};
  /// dz over the transverse length; dz/dr for the parabola.
  float cotTheta{0};
  /// Signed transverse curvature, positive clockwise seen from +z: a positive
  /// charge in a positive field moving along the fit direction.
  float curvature{0};
  /// q/p if the model knows the field, for a charge moving along the fit
  /// direction.
  std::optional<float> qOverP;

  /// The same track traversed the other way: a fit run in the opposite
  /// direction gives the reversed parameters.
  /// @return Parameters with direction and charge inverted
  SpacePointFilterParameters reversed() const {
    SpacePointFilterParameters out = *this;
    out.phi = Acts::detail::radian_sym(phi + std::numbers::pi_v<float>);
    out.cotTheta = -cotTheta;
    out.curvature = -curvature;
    if (qOverP.has_value()) {
      out.qOverP = -*qOverP;
    }
    return out;
  }
};

/// Kalman filter on space points for seeding.
///
/// `initialize` starts at one space point, in a frame towards a second one.
/// `predict` goes to the next space point and gives the residuals and their
/// chi2. `apply` takes the space point, and `parameters` gives the track at
/// the last one. Gates, cuts and scores belong to the caller, which decides
/// between `predict` and `apply`.
///
/// A filter family owns its state and its update. Inside a family the
/// propagation and the resolution are policies.
template <typename F>
concept SpacePointFilter = requires(const F& f, typename F::State& state,
                                    const SpacePointFilterInput& input,
                                    const typename F::Prediction& prediction) {
  { f.initialize(input, input) } -> std::same_as<typename F::State>;
  {
    f.predict(std::as_const(state), input, input)
  } -> std::same_as<typename F::Prediction>;
  { prediction.chi2Transverse() } -> std::same_as<float>;
  { prediction.chi2Longitudinal() } -> std::same_as<float>;
  f.apply(state, prediction);
  {
    f.parameters(std::as_const(state), input)
  } -> std::same_as<SpacePointFilterParameters>;
};

/// State of the decoupled family.
///
/// Frame: x' along the line from the first to the second space point, fixed
/// for the fit. Transverse state y', dy'/dx' and a bending variable, the
/// propagation says which. Longitudinal state z and dz/dl, the propagation
/// says what l is.
struct DecoupledSpacePointFilterState final {
  /// y', dy'/dx', bending variable
  std::array<float, 3> x{};
  /// z, dz/dl
  std::array<float, 2> y{};
  /// Upper triangle of the transverse covariance: 00 01 02 11 12 22
  std::array<float, 6> cx{};
  /// Upper triangle of the longitudinal covariance: 00 01 11
  std::array<float, 3> cy{};
  /// x' of the last space point
  float refX{};
  /// l of the last space point
  float refL{};
  /// Cosine of the frame
  float c{};
  /// Sine of the frame
  float s{};
};

/// One step of a decoupled propagation.
struct DecoupledSpacePointFilterStep final {
  /// x' step
  float a{};
  /// Coefficient of the bending variable in y'
  float b{};
  /// Coefficient of the bending variable in dy'/dx'
  float g{};
  /// l step
  float dl{};
  /// l at the measured space point
  float l{};
  /// Scattering variance of dy'/dx', added at the previous space point
  float varianceTransverse{};
  /// Scattering variance of dz/dl, added at the previous space point
  float varianceLongitudinal{};
};

namespace detail {

/// Transverse arc length for a chord on a circle of curvature `kappa`, to
/// second order in the turning angle.
inline float spacePointFilterArcLength(float chord, float kappa) {
  const float k2c2 = kappa * kappa * chord * chord;
  return chord * (1 + k2c2 / 24);
}

/// Highland scattering angle squared (PDG).
inline float spacePointFilterTheta0Squared(float pathInX0, float p,
                                           const ParticleHypothesis& particle) {
  assert(p > 0 && "momentum must be positive");
  // no material, or no bending measured yet
  if (pathInX0 <= 0 || std::isinf(p)) {
    return 0;
  }
  const float q = particle.absoluteCharge();
  const float beta = p / fastHypot(p, particle.mass());
  const float scale = 13.6 * UnitConstants::MeV;
  const float theta0 =
      scale * q / (beta * p) * std::sqrt(pathInX0) *
      (1 + 0.038f * std::log(pathInX0 * q * q / (beta * beta)));
  return theta0 * theta0;
}

/// Material crossed at a space point: a cylinder at 1/sin(theta), a disc at
/// 1/cos(theta).
inline float spacePointFilterPathInX0(const SpacePointFilterInput& input,
                                      float tau) {
  const float invSin = fastHypot(1.0f, tau);
  if (input.measuresR) {
    constexpr float minTau = 1e-3f;
    return input.thicknessInX0 * invSin / std::max(std::abs(tau), minTau);
  }
  return input.thicknessInX0 * invSin;
}

/// Cosine between the radial direction at the transverse position (x, y) and
/// the transverse direction of the fit.
inline float spacePointFilterRadialAlong(float x, float y, float c, float s,
                                         float slope) {
  const float tx = c - s * slope;
  const float ty = s + c * slope;
  return (x * tx + y * ty) / (fastHypot(tx, ty) * fastHypot(x, y));
}

/// Field bending the track in the transverse plane, Bz - tau * Br r . t: a
/// radial field along the track bends it through its motion along z.
inline float spacePointFilterBendingField(float br, float bz, float x, float y,
                                          float c, float s, float slope,
                                          float tau) {
  return bz - tau * br * spacePointFilterRadialAlong(x, y, c, s, slope);
}

/// Transverse distance between two space points.
inline float spacePointFilterChord(const SpacePointFilterInput& a,
                                   const SpacePointFilterInput& b) {
  return fastHypot(b.x - a.x, b.y - a.y);
}

/// Cosine and sine of the frame of a fit from `from` towards `to`.
inline std::array<float, 2> spacePointFilterFrame(
    const SpacePointFilterInput& from, const SpacePointFilterInput& to) {
  const float dx = to.x - from.x;
  const float dy = to.y - from.y;
  // not fastHypot: the GBTS tracking filter's dx * dx + dy * dy contracts to
  // an FMA in one expression, across fastHypot's calls it does not
  const float L = std::sqrt(dx * dx + dy * dy);
  assert(L > 0 && "space points must be apart in the transverse plane");
  return {dx / L, dy / L};
}

/// Position, direction and slope at the last space point.
inline SpacePointFilterParameters spacePointFilterGeometry(float refX, float c,
                                                           float s, float y,
                                                           float slope, float z,
                                                           float tau) {
  SpacePointFilterParameters out;
  out.position = Eigen::Vector3f(refX * c - y * s, refX * s + y * c, z);
  out.phi = std::atan2(s + c * slope, c - s * slope);
  out.cotTheta = tau;
  return out;
}

}  // namespace detail

/// Resolution from configured values: pixel resolutions, and for strips the
/// strip extent from the half vector.
///
/// Resolutions give 1 / (predicted variance + measurement variance) for the
/// transverse and the longitudinal residual. They take the predicted
/// variances so that the parabola filter keeps the arithmetic of the GBTS
/// tracking filter exactly: the compiler contracts `a * b + c` within one
/// expression.
struct ConfiguredResolution final {
  /// Configuration of the resolution.
  struct Config {
    /// Pixel resolution across the track.
    float sigmaX = 0.08 * UnitConstants::mm;
    /// Pixel resolution of the longitudinal coordinate.
    float sigmaY = 0.25 * UnitConstants::mm;
    /// Strip resolution across the strips.
    float sigmaXStrip = 0.03 * UnitConstants::mm;
  };

  /// Configuration
  Config cfg{};

  /// @param measured Space point of the measurement
  /// @param c Cosine of the frame
  /// @param s Sine of the frame
  /// @param tau Longitudinal slope
  /// @param predictedTransverse Predicted variance of the transverse position
  /// @param predictedLongitudinal Predicted variance of z
  /// @return Inverse residual variances, transverse and longitudinal
  std::array<float, 2> inverseResidualVariances(
      const SpacePointFilterInput& measured, const float c, const float s,
      const float tau, const float predictedTransverse,
      const float predictedLongitudinal) const {
    // expression for expression from GbtsTrackingFilter::update
    float sigmaX = cfg.sigmaX;
    float sigmaY = cfg.sigmaY;
    if (measured.stripHalfVector.has_value()) {
      constexpr float invSqrt3 = std::numbers::inv_sqrt3_v<float>;
      const std::array<float, 3>& half = *measured.stripHalfVector;
      const float alongX = -half[0] * s + half[1] * c;
      sigmaX = fastHypot(cfg.sigmaXStrip, alongX * invSqrt3);
      sigmaY = (measured.measuresR ? fastHypot(half[0], half[1])
                                   : std::abs(half[2])) *
               invSqrt3;
    }

    // a measured r becomes a z variance through the slope
    const float sigmaRZ =
        measured.measuresR ? square(sigmaY * tau) : sigmaY * sigmaY;

    return {1.0f / (predictedTransverse + sigmaX * sigmaX),
            1.0f / (predictedLongitudinal + sigmaRZ)};
  }
};

/// Resolution from the space point's own variances. A measured r becomes a z
/// variance through the longitudinal slope.
struct SpacePointResolution final {
  /// @copydoc ConfiguredResolution::inverseResidualVariances
  std::array<float, 2> inverseResidualVariances(
      const SpacePointFilterInput& measured, const float /*c*/,
      const float /*s*/, const float tau, const float predictedTransverse,
      const float predictedLongitudinal) const {
    const float varianceZ = measured.measuresR
                                ? measured.varianceLongitudinal * tau * tau
                                : measured.varianceLongitudinal;
    return {1.0f / (predictedTransverse + measured.varianceTransverse),
            1.0f / (predictedLongitudinal + varianceZ)};
  }
};

/// Parabola in the frame and a line in r-z, the GBTS tracking filter.
///
/// Bending variable d2y'/dx'2, l = r. Scattering from a fixed radiation
/// length, with the curvature as a momentum proxy.
struct ParabolaPropagation final {
  /// Configuration of the propagation.
  struct Config {
    /// Scattering angle at the reference curvature (900 MeV at eta = 0).
    float sigmaMS = 0.016;
    /// Radiation length per layer.
    float radLen = 0.025;
    /// Curvature at which the scattering angle is `sigmaMS`.
    float referenceCurvature = 1e-3 / UnitConstants::mm;
  };

  /// Configuration
  Config cfg{};

  /// @param from First space point
  /// @return l at the first space point
  float initialReference(const SpacePointFilterInput& from) const {
    return from.r;
  }

  /// @param from First space point
  /// @param to Second space point
  /// @return dz/dl from the two space points
  float initialLongitudinalSlope(const SpacePointFilterInput& from,
                                 const SpacePointFilterInput& to) const {
    return (to.z - from.z) / (to.r - from.r);
  }

  /// @param state State at the previous space point
  /// @param previous Previous space point
  /// @param measured Measured space point
  /// @param a x' step
  /// @return The step
  DecoupledSpacePointFilterStep step(
      const DecoupledSpacePointFilterState& state,
      const SpacePointFilterInput& previous,
      const SpacePointFilterInput& measured, const float a) const {
    // expression for expression from GbtsTrackingFilter::update
    const float tau2 = state.y[1] * state.y[1];
    const float invSin2 = 1 + tau2;
    const float lenCorr = previous.measuresR ? invSin2 / tau2 : invSin2;
    const float minPtFrac = std::abs(state.x[2]) / cfg.referenceCurvature;
    const float corrMS = cfg.sigmaMS * minPtFrac;
    const float sigma2 = cfg.radLen * lenCorr * corrMS * corrMS;

    DecoupledSpacePointFilterStep out;
    out.a = a;
    // in double, as the GBTS tracking filter evaluates it
    out.b = static_cast<float>(0.5 * a * a);
    out.g = a;
    out.l = measured.r;
    out.dl = out.l - state.refL;
    out.varianceTransverse = sigma2;
    out.varianceLongitudinal = sigma2;
    return out;
  }

  /// @param state State at the last space point
  /// @param last Last space point
  /// @return Signed transverse curvature
  float curvature(const DecoupledSpacePointFilterState& state,
                  const SpacePointFilterInput& /*last*/) const {
    const float slope2 = 1 + state.x[1] * state.x[1];
    return -state.x[2] / (slope2 * std::sqrt(slope2));
  }

  /// @param state State at the last space point
  /// @param last Last space point
  /// @return q/p, unknown without a field
  std::optional<float> qOverP(const DecoupledSpacePointFilterState& /*state*/,
                              const SpacePointFilterInput& /*last*/) const {
    return std::nullopt;
  }
};

/// Helix with the field per step.
///
/// Bending variable q/pT with Bz at the mean of the two space points, l =
/// transverse arc length. Highland scattering with p from the state and the
/// space point's material.
struct HelixPropagation final {
  /// Configuration of the propagation.
  struct Config {
    /// Particle hypothesis for the scattering: mass and absolute charge.
    ParticleHypothesis particle = ParticleHypothesis::pion();
  };

  /// Configuration
  Config cfg{};

  /// @copydoc ParabolaPropagation::initialReference
  float initialReference(const SpacePointFilterInput& /*from*/) const {
    return 0;
  }

  /// @copydoc ParabolaPropagation::initialLongitudinalSlope
  float initialLongitudinalSlope(const SpacePointFilterInput& from,
                                 const SpacePointFilterInput& to) const {
    return (to.z - from.z) / detail::spacePointFilterChord(from, to);
  }

  /// @copydoc ParabolaPropagation::step
  DecoupledSpacePointFilterStep step(
      const DecoupledSpacePointFilterState& state,
      const SpacePointFilterInput& previous,
      const SpacePointFilterInput& measured, const float a) const {
    const float qOverPt = state.x[2];
    const float slope = state.x[1];
    const float tau = state.y[1];
    // bending field at the middle of the step
    const float bz = detail::spacePointFilterBendingField(
        (previous.br + measured.br) / 2, (previous.bz + measured.bz) / 2,
        (previous.x + measured.x) / 2, (previous.y + measured.y) / 2, state.c,
        state.s, slope, tau);

    // y'' = q/pT * g, with g at the slope halfway through the step
    const float slope2 = 1 + slope * slope;
    const float g0 = -bz * slope2 * std::sqrt(slope2);
    const float mid = slope + qOverPt * g0 * a / 2;
    const float mid2 = 1 + mid * mid;
    const float g = -bz * mid2 * std::sqrt(mid2);

    // scattering at the space point we leave
    const float invSin2 = 1 + tau * tau;
    const float p =
        cfg.particle.absoluteCharge() * std::sqrt(invSin2) / std::abs(qOverPt);
    const float theta02 = detail::spacePointFilterTheta0Squared(
        detail::spacePointFilterPathInX0(previous, tau), p, cfg.particle);

    DecoupledSpacePointFilterStep out;
    out.a = a;
    out.b = g * a * a / 2;
    out.g = g * a;
    out.dl = detail::spacePointFilterArcLength(
        detail::spacePointFilterChord(previous, measured), qOverPt * bz);
    out.l = state.refL + out.dl;
    // d(dy'/dx') = (1 + y'^2) dphi, dphi = theta0 / sin(theta)
    out.varianceTransverse = theta02 * invSin2 * slope2 * slope2;
    // d(cot theta) = (1 + cot^2) dtheta
    out.varianceLongitudinal = theta02 * invSin2 * invSin2;
    return out;
  }

  /// @copydoc ParabolaPropagation::curvature
  float curvature(const DecoupledSpacePointFilterState& state,
                  const SpacePointFilterInput& last) const {
    return state.x[2] * detail::spacePointFilterBendingField(
                            last.br, last.bz, last.x, last.y, state.c, state.s,
                            state.x[1], state.y[1]);
  }

  /// @param state State at the last space point
  /// @param last Last space point
  /// @return q/p
  std::optional<float> qOverP(const DecoupledSpacePointFilterState& state,
                              const SpacePointFilterInput& /*last*/) const {
    const float tau = state.y[1];
    return state.x[2] / fastHypot(1.0f, tau);
  }
};

/// Prediction of the decoupled family.
struct DecoupledSpacePointFilterPrediction final {
  /// Predicted state at the measured space point
  DecoupledSpacePointFilterState state;
  /// Transverse residual
  float residualTransverse{};
  /// Longitudinal residual
  float residualLongitudinal{};
  /// Inverse variance of the transverse residual
  float inverseVarianceTransverse{};
  /// Inverse variance of the longitudinal residual
  float inverseVarianceLongitudinal{};

  /// @return chi2 of the transverse residual
  float chi2Transverse() const {
    return residualTransverse * residualTransverse * inverseVarianceTransverse;
  }
  /// @return chi2 of the longitudinal residual
  float chi2Longitudinal() const {
    return residualLongitudinal * residualLongitudinal *
           inverseVarianceLongitudinal;
  }
};

/// Decoupled family: 3 + 2 state, one-number updates, no inverse.
template <typename Propagation, typename Resolution>
class DecoupledSpacePointFilter final {
 public:
  /// State type
  using State = DecoupledSpacePointFilterState;
  /// Prediction type
  using Prediction = DecoupledSpacePointFilterPrediction;

  /// Configuration of the filter.
  struct Config {
    /// Initial variance of y', dy'/dx' and the bending variable.
    std::array<float, 3> initialVarianceTransverse{};
    /// Initial variance of z and dz/dl.
    std::array<float, 2> initialVarianceLongitudinal{};
  };

  /// @param config Configuration
  /// @param propagation Propagation policy
  /// @param resolution Resolution policy
  DecoupledSpacePointFilter(const Config& config,
                            const Propagation& propagation,
                            const Resolution& resolution)
      : m_cfg(config), m_propagation(propagation), m_resolution(resolution) {}

  /// State at `from`, in the frame along the line to `to`. The fit runs from
  /// `from` towards `to`.
  /// @param from First space point
  /// @param to Second space point
  /// @return Initial state
  State initialize(const SpacePointFilterInput& from,
                   const SpacePointFilterInput& to) const {
    // expression for expression from GbtsEdgeState::initialize
    State state;
    const auto [c, s] = detail::spacePointFilterFrame(from, to);
    state.c = c;
    state.s = s;

    state.refL = m_propagation.initialReference(from);
    state.refX = from.x * state.c + from.y * state.s;

    state.x[0] = -from.x * state.s + from.y * state.c;
    state.x[1] = 0;
    state.x[2] = 0;

    state.y[0] = from.z;
    state.y[1] = m_propagation.initialLongitudinalSlope(from, to);

    state.cx[0] = m_cfg.initialVarianceTransverse[0];
    state.cx[3] = m_cfg.initialVarianceTransverse[1];
    state.cx[5] = m_cfg.initialVarianceTransverse[2];

    state.cy[0] = m_cfg.initialVarianceLongitudinal[0];
    state.cy[2] = m_cfg.initialVarianceLongitudinal[1];
    return state;
  }

  /// Predict to `measured`, scattering at `previous`.
  /// @param state State at the previous space point
  /// @param previous Previous space point
  /// @param measured Measured space point
  /// @return Prediction with the residuals
  Prediction predict(const State& state, const SpacePointFilterInput& previous,
                     const SpacePointFilterInput& measured) const {
    // expression for expression from GbtsTrackingFilter::update; the
    // coefficient of the bending variable in dy'/dx' is `G`, `A` there
    const float refX = measured.x * state.c + measured.y * state.s;
    const float mx = -measured.x * state.s + measured.y * state.c;

    const DecoupledSpacePointFilterStep st =
        m_propagation.step(state, previous, measured, refX - state.refX);

    const float cx00 = state.cx[0];
    const float cx01 = state.cx[1];
    const float cx02 = state.cx[2];
    const float cx11 = state.cx[3] + st.varianceTransverse;
    const float cx12 = state.cx[4];
    const float cx22 = state.cx[5];
    const float cy00 = state.cy[0];
    const float cy01 = state.cy[1];
    const float cy11 = state.cy[2] + st.varianceLongitudinal;

    const float A = st.a;
    const float B = st.b;
    const float G = st.g;
    const float dl = st.dl;

    Prediction p;
    State& out = p.state;
    out.c = state.c;
    out.s = state.s;
    out.refX = refX;
    out.refL = st.l;

    out.x[0] = state.x[0] + state.x[1] * A + state.x[2] * B;
    out.x[1] = state.x[1] + state.x[2] * G;
    out.x[2] = state.x[2];

    out.cx[0] = cx00 + 2 * cx01 * A + 2 * cx02 * B + A * A * cx11 +
                2 * A * B * cx12 + B * B * cx22;
    out.cx[1] =
        cx01 + cx11 * A + cx12 * B + cx02 * G + A * G * cx12 + G * B * cx22;
    out.cx[2] = cx02 + cx12 * A + cx22 * B;
    out.cx[3] = cx11 + 2 * G * cx12 + G * G * cx22;
    out.cx[4] = cx12 + cx22 * G;
    out.cx[5] = cx22;

    out.y[0] = state.y[0] + state.y[1] * dl;
    out.y[1] = state.y[1];

    out.cy[0] = cy00 + 2 * cy01 * dl + dl * dl * cy11;
    out.cy[1] = cy01 + dl * cy11;
    out.cy[2] = cy11;

    p.residualTransverse = mx - out.x[0];
    p.residualLongitudinal = measured.z - out.y[0];

    const auto [dx, dy] = m_resolution.inverseResidualVariances(
        measured, state.c, state.s, state.y[1], out.cx[0], out.cy[0]);
    p.inverseVarianceTransverse = dx;
    p.inverseVarianceLongitudinal = dy;
    return p;
  }

  /// Take the measured space point of a prediction.
  /// @param state State to update
  /// @param prediction Prediction to the measured space point
  void apply(State& state, const Prediction& prediction) const {
    const State& pr = prediction.state;
    const float Dx = prediction.inverseVarianceTransverse;
    const float Dy = prediction.inverseVarianceLongitudinal;

    const std::array<float, 3> CHx = {pr.cx[0], pr.cx[1], pr.cx[2]};
    const std::array<float, 2> CHy = {pr.cy[0], pr.cy[1]};
    const std::array<float, 3> Kx = {Dx * CHx[0], Dx * CHx[1], Dx * CHx[2]};
    const std::array<float, 2> Ky = {Dy * CHy[0], Dy * CHy[1]};

    state.c = pr.c;
    state.s = pr.s;
    state.refX = pr.refX;
    state.refL = pr.refL;

    for (std::size_t i = 0; i < 3; ++i) {
      state.x[i] = pr.x[i] + Kx[i] * prediction.residualTransverse;
    }
    for (std::size_t i = 0; i < 2; ++i) {
      state.y[i] = pr.y[i] + Ky[i] * prediction.residualLongitudinal;
    }

    // upper triangles, entry (i, j) as Cx[i][j] - Kx[i] * CHx[j]
    state.cx[0] = pr.cx[0] - Kx[0] * CHx[0];
    state.cx[1] = pr.cx[1] - Kx[0] * CHx[1];
    state.cx[2] = pr.cx[2] - Kx[0] * CHx[2];
    state.cx[3] = pr.cx[3] - Kx[1] * CHx[1];
    state.cx[4] = pr.cx[4] - Kx[1] * CHx[2];
    state.cx[5] = pr.cx[5] - Kx[2] * CHx[2];

    state.cy[0] = pr.cy[0] - Ky[0] * CHy[0];
    state.cy[1] = pr.cy[1] - Ky[0] * CHy[1];
    state.cy[2] = pr.cy[2] - Ky[1] * CHy[1];
  }

  /// @param state State at the last space point
  /// @param last Last space point
  /// @return Track parameters at the last space point
  SpacePointFilterParameters parameters(
      const State& state, const SpacePointFilterInput& last) const {
    SpacePointFilterParameters out = detail::spacePointFilterGeometry(
        state.refX, state.c, state.s, state.x[0], state.x[1], state.y[0],
        state.y[1]);
    out.curvature = m_propagation.curvature(state, last);
    out.qOverP = m_propagation.qOverP(state, last);
    return out;
  }

  /// @return Configuration
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  Propagation m_propagation;
  Resolution m_resolution;
};

/// State of the coupled family: y', dy'/dx', z, dz/ds, q/p in the decoupled
/// family's frame, l = transverse arc length.
struct CoupledSpacePointFilterState final {
  /// y', dy'/dx', z, dz/ds, q/p
  Eigen::Matrix<float, 5, 1> x{Eigen::Matrix<float, 5, 1>::Zero()};
  /// Covariance
  Eigen::Matrix<float, 5, 5> cov{Eigen::Matrix<float, 5, 5>::Zero()};
  /// x' of the last space point
  float refX{};
  /// l of the last space point
  float refL{};
  /// Cosine of the frame
  float c{};
  /// Sine of the frame
  float s{};
};

/// Prediction of the coupled family.
struct CoupledSpacePointFilterPrediction final {
  /// Predicted state at the measured space point
  CoupledSpacePointFilterState state;
  /// The measured space point, for the second residual after the first update
  SpacePointFilterInput measured;
  /// Transverse residual
  float residualTransverse{};
  /// Longitudinal residual
  float residualLongitudinal{};
  /// Inverse variance of the transverse residual
  float inverseVarianceTransverse{};
  /// Inverse variance of the longitudinal residual
  float inverseVarianceLongitudinal{};

  /// @return Marginal chi2 of the transverse residual, before either update
  float chi2Transverse() const {
    return residualTransverse * residualTransverse * inverseVarianceTransverse;
  }
  /// @return Marginal chi2 of the longitudinal residual, before either update
  float chi2Longitudinal() const {
    return residualLongitudinal * residualLongitudinal *
           inverseVarianceLongitudinal;
  }
};

/// Coupled family: helix with q/p.
///
/// The bending depends on the longitudinal slope, so the planes share a
/// covariance. The two residuals are taken one after the other, which is
/// exact for independent measurement errors and needs no inverse.
template <typename Resolution>
class CoupledSpacePointFilter final {
 public:
  /// State type
  using State = CoupledSpacePointFilterState;
  /// Prediction type
  using Prediction = CoupledSpacePointFilterPrediction;

  /// Configuration of the filter.
  struct Config {
    /// Initial variance of y', dy'/dx', z, dz/ds, q/p.
    std::array<float, 5> initialVariance{};
    /// Particle hypothesis for the scattering: mass and absolute charge.
    ParticleHypothesis particle = ParticleHypothesis::pion();
  };

  /// @param config Configuration
  /// @param resolution Resolution policy
  CoupledSpacePointFilter(const Config& config, const Resolution& resolution)
      : m_cfg(config), m_resolution(resolution) {}

  /// @copydoc DecoupledSpacePointFilter::initialize
  State initialize(const SpacePointFilterInput& from,
                   const SpacePointFilterInput& to) const {
    State state;
    const auto [c, s] = detail::spacePointFilterFrame(from, to);
    state.c = c;
    state.s = s;
    state.refX = from.x * c + from.y * s;

    state.x(0) = -from.x * s + from.y * c;
    state.x(1) = 0;
    state.x(2) = from.z;
    state.x(3) = (to.z - from.z) / detail::spacePointFilterChord(from, to);
    state.x(4) = 0;

    for (int i = 0; i < 5; ++i) {
      state.cov(i, i) = m_cfg.initialVariance[i];
    }
    return state;
  }

  /// @copydoc DecoupledSpacePointFilter::predict
  Prediction predict(const State& state, const SpacePointFilterInput& previous,
                     const SpacePointFilterInput& measured) const {
    using Matrix5 = Eigen::Matrix<float, 5, 5>;

    const float refX = measured.x * state.c + measured.y * state.s;
    const float mx = -measured.x * state.s + measured.y * state.c;
    const float a = refX - state.refX;

    const float slope = state.x(1);
    const float tau = state.x(3);
    const float qOverP = state.x(4);
    // bending field at the middle of the step, Bz - tau * bAlong
    const float bAlong =
        (previous.br + measured.br) / 2 *
        detail::spacePointFilterRadialAlong((previous.x + measured.x) / 2,
                                            (previous.y + measured.y) / 2,
                                            state.c, state.s, slope);
    const float bz = (previous.bz + measured.bz) / 2 - tau * bAlong;

    // y'' = q/p * g with g = -B / sin(theta) * (1 + y'^2)^(3/2), at the slope
    // halfway through the step; B depends on tau, its dependence on the slope
    // through the direction is left out
    const float invSin2 = 1 + tau * tau;
    const float invSin = std::sqrt(invSin2);
    const float slope2 = 1 + slope * slope;
    const float g0 = -bz * invSin * slope2 * std::sqrt(slope2);
    const float mid = slope + qOverP * g0 * a / 2;
    const float mid2 = 1 + mid * mid;
    const float midRoot = std::sqrt(mid2);
    const float g = -bz * invSin * mid2 * midRoot;
    const float dgdSlope = -bz * invSin * 3 * mid * midRoot;
    const float dgdTau = (bAlong * invSin - bz * tau / invSin) * mid2 * midRoot;

    const float ds = detail::spacePointFilterArcLength(
        detail::spacePointFilterChord(previous, measured),
        qOverP * invSin * bz);
    const float a2 = a * a / 2;

    // scattering at the space point we leave
    const float theta02 = detail::spacePointFilterTheta0Squared(
        detail::spacePointFilterPathInX0(previous, tau),
        m_cfg.particle.absoluteCharge() / std::abs(qOverP), m_cfg.particle);
    Matrix5 cov = state.cov;
    cov(1, 1) += theta02 * invSin2 * slope2 * slope2;
    cov(3, 3) += theta02 * invSin2 * invSin2;

    Matrix5 F = Matrix5::Identity();
    F(0, 1) = a + qOverP * a2 * dgdSlope;
    F(0, 3) = qOverP * a2 * dgdTau;
    F(0, 4) = g * a2;
    F(1, 1) = 1 + qOverP * a * dgdSlope;
    F(1, 3) = qOverP * a * dgdTau;
    F(1, 4) = g * a;
    F(2, 3) = ds;

    Prediction p;
    State& out = p.state;
    out.c = state.c;
    out.s = state.s;
    out.refX = refX;
    out.refL = state.refL + ds;

    out.x = state.x;
    out.x(0) += slope * a + qOverP * g * a2;
    out.x(1) += qOverP * g * a;
    out.x(2) += tau * ds;
    out.cov = F * cov * F.transpose();

    p.measured = measured;
    p.residualTransverse = mx - out.x(0);
    p.residualLongitudinal = measured.z - out.x(2);

    const auto [dx, dy] = m_resolution.inverseResidualVariances(
        measured, state.c, state.s, tau, out.cov(0, 0), out.cov(2, 2));
    p.inverseVarianceTransverse = dx;
    p.inverseVarianceLongitudinal = dy;
    return p;
  }

  /// @copydoc DecoupledSpacePointFilter::apply
  void apply(State& state, const Prediction& prediction) const {
    using Vector5 = Eigen::Matrix<float, 5, 1>;

    state = prediction.state;

    // transverse residual
    {
      const Vector5 k = state.cov.col(0) * prediction.inverseVarianceTransverse;
      state.x += k * prediction.residualTransverse;
      state.cov -= k * state.cov.row(0);
    }

    // longitudinal residual, after the transverse update
    {
      const float residual = prediction.measured.z - state.x(2);
      const float d = m_resolution.inverseResidualVariances(
          prediction.measured, state.c, state.s, state.x(3), state.cov(0, 0),
          state.cov(2, 2))[1];
      const Vector5 k = state.cov.col(2) * d;
      state.x += k * residual;
      state.cov -= k * state.cov.row(2);
    }

    state.cov = (state.cov + state.cov.transpose()) / 2;
  }

  /// @copydoc DecoupledSpacePointFilter::parameters
  SpacePointFilterParameters parameters(
      const State& state, const SpacePointFilterInput& last) const {
    SpacePointFilterParameters out = detail::spacePointFilterGeometry(
        state.refX, state.c, state.s, state.x(0), state.x(1), state.x(2),
        state.x(3));
    const float qOverP = state.x(4);
    out.curvature = qOverP * fastHypot(1.0f, out.cotTheta) *
                    detail::spacePointFilterBendingField(
                        last.br, last.bz, last.x, last.y, state.c, state.s,
                        state.x(1), state.x(3));
    out.qOverP = qOverP;
    return out;
  }

  /// @return Configuration
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  Resolution m_resolution;
};

/// Parabola filter, the GBTS tracking filter.
using ParabolaSpacePointFilter =
    DecoupledSpacePointFilter<ParabolaPropagation, ConfiguredResolution>;
/// Decoupled helix with q/pT, the field per step and Highland scattering.
using HelixSpacePointFilter =
    DecoupledSpacePointFilter<HelixPropagation, SpacePointResolution>;
/// Coupled helix with q/p.
using CoupledHelixSpacePointFilter =
    CoupledSpacePointFilter<SpacePointResolution>;

static_assert(SpacePointFilter<ParabolaSpacePointFilter>);
static_assert(SpacePointFilter<HelixSpacePointFilter>);
static_assert(SpacePointFilter<CoupledHelixSpacePointFilter>);

}  // namespace Acts::Experimental
