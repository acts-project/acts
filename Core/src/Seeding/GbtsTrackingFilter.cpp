// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsTrackingFilter.hpp"

#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <utility>

namespace Acts::Experimental::detail {

void detail::GbtsEdgeState::initialize(const detail::GbtsEdge& pS,
                                       const detail::GbtsNodeView& nodeView,
                                       const std::array<float, 3>& varianceX,
                                       const std::array<float, 2>& varianceY) {
  j = 0;

  // The chain is followed inwards, from the outer node n2 of the head edge to
  // its inner node n1.
  const detail::GbtsNodeProxy n1 = nodeView[pS.n1];
  const detail::GbtsNodeProxy n2 = nodeView[pS.n2];

  const float dx = n1.x() - n2.x();
  const float dy = n1.y() - n2.y();
  const float L = std::sqrt(dx * dx + dy * dy);

  // The transverse state lives in a frame rotated so that x' points along the
  // head edge:
  //   x' =  x c + y s
  //   y' = -x s + y c
  s = dy / L;
  c = dx / L;

  refY = n2.r();
  refX = n2.x() * c + n2.y() * s;

  // Transverse state at n2: y', the slope dy'/dx' (zero along the head edge)
  // and the curvature (not known yet).
  x[0] = -n2.x() * s + n2.y() * c;
  x[1] = 0;
  x[2] = 0;

  // Longitudinal state at n2: z and the slope dz/ds, with s the transverse
  // arc length, growing outwards. Approximation: the slope is taken over the
  // chord of the head edge, as the curvature is not known yet.
  y[0] = n2.z();
  y[1] = (n1.z() - n2.z()) / std::copysign(L, n1.r() - n2.r());

  cx = {varianceX[0], 0, 0, varianceX[1], 0, varianceX[2]};

  cy = {varianceY[0], 0, varianceY[1]};
}

}  // namespace Acts::Experimental::detail

namespace Acts::Experimental {

GbtsTrackingFilter::GbtsTrackingFilter(
    const Config& config, const std::shared_ptr<const GbtsGeometry>& geometry,
    std::unique_ptr<const Logger> logger)
    : m_cfg(config), m_geometry(geometry), m_logger(std::move(logger)) {
  const float highland = 0.0136f * Acts::UnitConstants::GeV *
                         (1 + 0.038f * std::log(m_cfg.radLen));
  m_scattering = highland * highland * m_cfg.radLen;
  m_maxCurvature2 = m_cfg.maxCurvature * m_cfg.maxCurvature;
}

bool GbtsTrackingFilter::update(const detail::GbtsNodeView& nodeView,
                                const detail::GbtsEdge& pS,
                                detail::GbtsEdgeState& ts,
                                const float bFieldInZ) const {
  if (ts.cx[5] < 0 || ts.cx[3] < 0 || ts.cx[0] < 0) {
    ACTS_DEBUG("Negative cov_x");
  }

  if (ts.cy[2] < 0 || ts.cy[0] < 0) {
    ACTS_DEBUG("Negative cov_y");
  }

  // Multiple scattering in the layer of n2, the node the state is at, added
  // to the slope variances before the step. Approximations: the whole layer
  // acts at the node, the particle moves at the speed of light, and the
  // logarithmic term of the Highland formula uses the material at normal
  // incidence instead of the material along the path.

  const float tau2 = ts.y[1] * ts.y[1];
  const float invSin2 = 1 + tau2;

  const detail::GbtsNodeProxy n1 = nodeView[pS.n1];
  const detail::GbtsNodeProxy n2 = nodeView[pS.n2];

  const GbtsLayerType layerType1 = getLayerType(n2.layer());

  // Path through the layer over its thickness, squared: 1 / sin^2(theta) for
  // a barrel layer, 1 / cos^2(theta) for an endcap layer. Approximation: the
  // angle at which the track crosses the layer in the transverse plane is
  // neglected.
  const float lenCorr =
      layerType1 == GbtsLayerType::Barrel ? invSin2 : invSin2 / tau2;
  const float pathFactor = std::sqrt(lenCorr);

  // As long as the curvature is not known to better than maxCurvature, the
  // state scatters like the softest track that can pass the cut.
  float curvature = std::abs(ts.x[2]);
  if (ts.cx[5] > m_maxCurvature2) {
    curvature = std::max(curvature, m_cfg.maxCurvature);
  }

  // Variance of the azimuth from scattering: the Highland angle squared, with
  // 1/p^2 = curvature^2 / (B^2 (1 + tau^2)), times (1 + tau^2) because a
  // deflection by theta0 turns the azimuth by theta0 / sin(theta).
  const float slope2 = 1 + ts.x[1] * ts.x[1];
  const float theta2 = m_scattering * curvature * curvature * pathFactor /
                       (bFieldInZ * bFieldInZ);

  // The angles enter the slopes as d(dy'/dx') = (1 + (dy'/dx')^2) dphi and
  // d(dz/ds) = -(1 + tau^2) dtheta.
  ts.cx[3] += theta2 * slope2 * slope2;

  ts.cy[2] += theta2 * invSin2;

  // Prediction of the state at n1.

  std::array<float, 3> X{};
  std::array<float, 2> Y{};
  std::array<float, 6> Cx{};
  std::array<float, 3> Cy{};

  const float x = n1.x();
  const float y = n1.y();
  const float z = n1.z();
  const float r = n1.r();

  const float refX = x * ts.c + y * ts.s;
  const float mx = -x * ts.s + y * ts.c;  // measured y' of n1
  const float refY = r;
  const float my = z;  // measured z of n1

  // step along x'
  const float A = refX - ts.refX;

  // On a circle the sine of the direction phi against x' grows linearly with
  // x': sin(phi1) = sin(phi0) + curvature A.
  const float cosPhi0 = 1 / std::sqrt(slope2);
  const float sinPhi0 = ts.x[1] * cosPhi0;
  const float sinPhi1 = sinPhi0 + ts.x[2] * A;
  // A track turning more than 90 degrees against the head edge before it
  // reaches n1 has no point at the x' of n1, so the step is rejected.
  if (std::abs(sinPhi1) >= 1) {
    return false;
  }
  const float cosPhi1 = std::sqrt(1 - sinPhi1 * sinPhi1);
  const float invCosPhi1 = 1 / cosPhi1;
  // The change in y' is (cos(phi0) - cos(phi1)) / curvature, which equals A
  // times the tangent of the mean direction (phi0 + phi1) / 2. This form also
  // holds for a straight line.
  const float tanMean = (sinPhi0 + sinPhi1) / (cosPhi0 + cosPhi1);
  X[0] = ts.x[0] + A * tanMean;
  X[1] = sinPhi1 * invCosPhi1;
  X[2] = ts.x[2];

  // z advances with the transverse arc length between the nodes. The arc over
  // the chord, 2 asin(chord curvature / 2) / curvature, is expanded to second
  // order in chord * curvature. Its sign follows the step in r.
  const float chord = fastHypot(n1.x() - n2.x(), n1.y() - n2.y());
  const float arc =
      chord * (1 + chord * chord * ts.x[2] * ts.x[2] * (1.f / 24));
  const float dr = std::copysign(arc, refY - ts.refY);

  // Jacobian of y' and dy'/dx' after the step with respect to dy'/dx' and the
  // curvature before it, using dphi0/d(dy'/dx') = cos^2(phi0),
  // dphi1/d(dy'/dx') = cos^3(phi0) / cos(phi1) and
  // dphi1/dcurvature = A / cos(phi1).
  const float cos2Phi0 = cosPhi0 * cosPhi0;
  const float invCos3Phi1 = invCosPhi1 * invCosPhi1 * invCosPhi1;
  const float halfStepSec2 = 0.5f * A * (1 + tanMean * tanMean);
  const float j01 = halfStepSec2 * cos2Phi0 * (1 + cosPhi0 * invCosPhi1);
  const float j02 = halfStepSec2 * A * invCosPhi1;
  const float j11 = cos2Phi0 * cosPhi0 * invCos3Phi1;
  const float j12 = A * invCos3Phi1;

  // Predicted covariance from the Jacobian. The symmetric covariances are
  // stored as upper triangles: 00, 01, 02, 11, 12, 22.
  const auto& C = ts.cx;
  const float u00 = C[0] + j01 * C[1] + j02 * C[2];
  const float u10 = C[1] + j01 * C[3] + j02 * C[4];
  const float u20 = C[2] + j01 * C[4] + j02 * C[5];
  const float u01 = j11 * C[1] + j12 * C[2];
  const float u11 = j11 * C[3] + j12 * C[4];
  const float u21 = j11 * C[4] + j12 * C[5];
  Cx[0] = u00 + j01 * u10 + j02 * u20;
  Cx[1] = u01 + j01 * u11 + j02 * u21;
  Cx[2] = u20;
  Cx[3] = j11 * u11 + j12 * u21;
  Cx[4] = u21;
  Cx[5] = C[5];

  // Longitudinal prediction: z moves by the slope times the arc length.
  Y[0] = ts.y[0] + ts.y[1] * dr;
  Y[1] = ts.y[1];

  Cy[0] = ts.cy[0] + 2 * ts.cy[1] * dr + dr * dr * ts.cy[2];
  Cy[1] = ts.cy[1] + dr * ts.cy[2];
  Cy[2] = ts.cy[2];

  // Residuals and measurement errors.

  const float resid_x = mx - X[0];
  const float resid_y = my - Y[0];

  const std::array<float, 3> CHx = {Cx[0], Cx[1], Cx[2]};
  const std::array<float, 2> CHy = {Cy[0], Cy[1]};

  const GbtsLayerType type = getLayerType(n1.layer());

  // Strip error along the outer strip: uniform over its length, half length /
  // sqrt(3).
  float sigmaX = m_cfg.sigmaX;
  float sigmaY = m_cfg.sigmaY;
  if (const auto* strip = nodeView.strip(n1.index()); strip != nullptr) {
    constexpr float invSqrt3 = std::numbers::inv_sqrt3_v<float>;
    const std::array<float, 3>& half = strip->outerHalfVector;
    // the part of the strip length along y'
    const float alongX = -half[0] * ts.s + half[1] * ts.c;
    sigmaX = fastHypot(m_cfg.sigmaXStrip, alongX * invSqrt3);
    sigmaY = (type == GbtsLayerType::Barrel ? std::abs(half[2])
                                            : fastHypot(half[0], half[1])) *
             invSqrt3;
  }

  // The variance of the measured z is varY / normY (normY = 1 in the barrel),
  // kept as two factors so that the update below needs one division.
  float normY = 1;
  float varY = 0;
  if (type == GbtsLayerType::Barrel) {
    // a barrel node measures z
    varY = sigmaY * sigmaY;
  } else if (type == GbtsLayerType::Endcap) {
    // An endcap node has a known z and measures r. The error in r turns into
    // an error in z through dz/dr = (dz/ds) / cos(psi), with psi the angle
    // between the track direction at n1 and the radius.
    const float dirX = ts.c * cosPhi1 - ts.s * sinPhi1;
    const float dirY = ts.s * cosPhi1 + ts.c * sinPhi1;
    const float rCosPsi = x * dirX + y * dirY;
    const float sigmaZr = sigmaY * ts.y[1] * r;
    normY = rCosPsi * rCosPsi;
    varY = sigmaZr * sigmaZr;
  } else {
    throw std::runtime_error("invalid layer type");
  }

  const float varX = sigmaX * sigmaX;
  const float Dx = 1.0f / (Cx[0] + varX);

  const float invDenY = 1.0f / (normY * Cy[0] + varY);
  const float Dy = normY * invDenY;

  // A step failing either chi2 cut ends the branch.
  const float dchi2_x = resid_x * resid_x * Dx;
  const float dchi2_y = resid_y * resid_y * Dy;

  if (dchi2_x > m_cfg.maxDChi2X || dchi2_y > m_cfg.maxDChi2Y) {
    return false;
  }

  ts.j += m_cfg.addHit - dchi2_x * m_cfg.weightX - dchi2_y * m_cfg.weightY;

  // Kalman update of the state.

  const std::array<float, 3> Kx = {Dx * Cx[0], Dx * Cx[1], Dx * Cx[2]};
  const std::array<float, 2> Ky = {Dy * Cy[0], Dy * Cy[1]};

  for (std::uint32_t i = 0; i < 3; ++i) {
    ts.x[i] = X[i] + Kx[i] * resid_x;
  }

  if (std::abs(ts.x[2]) > m_cfg.maxCurvature) {
    return false;
  }

  for (std::uint32_t i = 0; i < 2; ++i) {
    ts.y[i] = Y[i] + Ky[i] * resid_y;
  }

  // z at the beam line. Approximation: extrapolated with the radius instead of
  // the arc length to the beam line.
  const float z0 = ts.y[0] - refY * ts.y[1];

  if (std::abs(z0) > m_cfg.maxZ0) {
    return false;
  }

  // Covariance after the update, C - K C. The measured row is written as
  // C R / (C00 + R) to avoid a subtraction that loses precision in float.
  const float keepX = varX * Dx;
  ts.cx[0] = Cx[0] * keepX;
  ts.cx[1] = Cx[1] * keepX;
  ts.cx[2] = Cx[2] * keepX;
  ts.cx[3] = Cx[3] - Kx[1] * CHx[1];
  ts.cx[4] = Cx[4] - Kx[1] * CHx[2];
  ts.cx[5] = Cx[5] - Kx[2] * CHx[2];

  const float keepY = varY * invDenY;
  ts.cy[0] = Cy[0] * keepY;
  ts.cy[1] = Cy[1] * keepY;
  ts.cy[2] = Cy[2] - Ky[1] * CHy[1];
  ts.refX = refX;
  ts.refY = refY;
  return true;
}

GbtsLayerType GbtsTrackingFilter::getLayerType(
    const GbtsLayerIndex layerIndex) const {
  return m_geometry->layerDescription(layerIndex).type;
}

}  // namespace Acts::Experimental
