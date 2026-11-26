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

  // n2->n1

  const detail::GbtsNodeProxy n1 = nodeView[pS.n1];
  const detail::GbtsNodeProxy n2 = nodeView[pS.n2];

  const float dx = n1.x() - n2.x();
  const float dy = n1.y() - n2.y();
  const float length = std::sqrt(dx * dx + dy * dy);

  s = dy / length;
  c = dx / length;

  // transform for extrapolation and update
  //  x' =  x*c + y*s
  //  y' = -x*s + y*c

  refY = n2.r();
  refX = n2.x() * c + n2.y() * s;

  // X-state: y, dy/dx, d2y/dx2

  x[0] = -n2.x() * s + n2.y() * c;
  x[1] = 0;
  x[2] = 0;

  // Y-state: z, dz/dr

  y[0] = n2.z();
  y[1] = (n1.z() - n2.z()) / (n1.r() - n2.r());

  cx = {};
  cx[0][0] = varianceX[0];
  cx[1][1] = varianceX[1];
  cx[2][2] = varianceX[2];

  cy = {};
  cy[0][0] = varianceY[0];
  cy[1][1] = varianceY[1];
}

}  // namespace Acts::Experimental::detail

namespace Acts::Experimental {

GbtsTrackingFilter::GbtsTrackingFilter(
    const Config& config, const std::shared_ptr<const GbtsGeometry>& geometry,
    std::unique_ptr<const Logger> logger)
    : m_cfg(config), m_geometry(geometry), m_logger(std::move(logger)) {}

bool GbtsTrackingFilter::update(const detail::GbtsNodeView& nodeView,
                                const detail::GbtsEdge& pS,
                                detail::GbtsEdgeState& ts) const {
  if (ts.cx[2][2] < 0 || ts.cx[1][1] < 0 || ts.cx[0][0] < 0) {
    ACTS_DEBUG("Negative cov_x");
  }

  if (ts.cy[1][1] < 0 || ts.cy[0][0] < 0) {
    ACTS_DEBUG("Negative cov_y");
  }

  // add ms.

  const float tau2 = ts.y[1] * ts.y[1];
  const float invSin2 = 1 + tau2;

  const detail::GbtsNodeProxy n1 = nodeView[pS.n1];
  const detail::GbtsNodeProxy n2 = nodeView[pS.n2];

  const GbtsLayerType layerType1 = getLayerType(n2.layer());

  const float lenCorr =
      layerType1 == GbtsLayerType::Barrel ? invSin2 : invSin2 / tau2;

  const float minPtFrac = std::abs(ts.x[2]) / m_cfg.maxCurvature;

  const float corrMS = m_cfg.sigmaMS * minPtFrac;

  const float sigma2 = m_cfg.radLen * lenCorr * corrMS * corrMS;  // /invSin2

  ts.cx[1][1] += sigma2;

  ts.cy[1][1] += sigma2;

  // extrapolation

  std::array<float, 3> predictedX{};
  std::array<float, 2> predictedY{};
  std::array<std::array<float, 3>, 3> predictedCovX{};
  std::array<std::array<float, 2>, 2> predictedCovY{};

  const float x = n1.x();
  const float y = n1.y();
  const float z = n1.z();
  const float r = n1.r();

  const float refX = x * ts.c + y * ts.s;
  const float mx = -x * ts.s + y * ts.c;  // measurement of the transverse state
  const float refY = r;
  const float my = z;  // measurement of the longitudinal state

  const float deltaX = refX - ts.refX;
  const float halfDeltaX2 = 0.5 * deltaX * deltaX;
  const float dr = refY - ts.refY;

  predictedX[0] = ts.x[0] + ts.x[1] * deltaX + ts.x[2] * halfDeltaX2;
  predictedX[1] = ts.x[1] + ts.x[2] * deltaX;
  predictedX[2] = ts.x[2];

  predictedCovX[0][0] =
      ts.cx[0][0] + 2 * ts.cx[0][1] * deltaX + 2 * ts.cx[0][2] * halfDeltaX2 +
      deltaX * deltaX * ts.cx[1][1] + 2 * deltaX * halfDeltaX2 * ts.cx[1][2] +
      halfDeltaX2 * halfDeltaX2 * ts.cx[2][2];
  predictedCovX[0][1] = predictedCovX[1][0] =
      ts.cx[0][1] + ts.cx[1][1] * deltaX + ts.cx[1][2] * halfDeltaX2 +
      ts.cx[0][2] * deltaX + deltaX * deltaX * ts.cx[1][2] +
      deltaX * halfDeltaX2 * ts.cx[2][2];
  predictedCovX[0][2] = predictedCovX[2][0] =
      ts.cx[0][2] + ts.cx[1][2] * deltaX + ts.cx[2][2] * halfDeltaX2;

  predictedCovX[1][1] =
      ts.cx[1][1] + 2 * deltaX * ts.cx[1][2] + deltaX * deltaX * ts.cx[2][2];
  predictedCovX[1][2] = predictedCovX[2][1] =
      ts.cx[1][2] + ts.cx[2][2] * deltaX;

  predictedCovX[2][2] = ts.cx[2][2];

  predictedY[0] = ts.y[0] + ts.y[1] * dr;
  predictedY[1] = ts.y[1];

  predictedCovY[0][0] =
      ts.cy[0][0] + 2 * ts.cy[0][1] * dr + dr * dr * ts.cy[1][1];
  predictedCovY[0][1] = predictedCovY[1][0] = ts.cy[0][1] + dr * ts.cy[1][1];
  predictedCovY[1][1] = ts.cy[1][1];

  // chi2 test

  const float residualX = mx - predictedX[0];
  const float residualY = my - predictedY[0];

  const std::array<float, 3> measurementCovX = {
      predictedCovX[0][0], predictedCovX[0][1], predictedCovX[0][2]};
  const std::array<float, 2> measurementCovY = {predictedCovY[0][0],
                                                predictedCovY[0][1]};

  float varianceRZ = 0;

  const GbtsLayerType type = getLayerType(n1.layer());

  // Across the strips a pair is sharper than a pixel, along them far worse:
  // the crossing is unconstrained over the strip, a half length over sqrt(3)
  // off its half vector -- the outer strip's, a stereo angle from the inner.
  float sigmaX = m_cfg.sigmaX;
  float sigmaY = m_cfg.sigmaY;
  if (const auto* strip = nodeView.strip(n1.index()); strip != nullptr) {
    constexpr float invSqrt3 = std::numbers::inv_sqrt3_v<float>;
    const std::array<float, 3>& half = strip->outerHalfVector;
    // the walk is along the strip, so project its reach onto the fit's axes
    const float alongX = -half[0] * ts.s + half[1] * ts.c;
    sigmaX = fastHypot(m_cfg.sigmaXStrip, alongX * invSqrt3);
    sigmaY = (type == GbtsLayerType::Barrel ? std::abs(half[2])
                                            : fastHypot(half[0], half[1])) *
             invSqrt3;
  }

  if (type == GbtsLayerType::Barrel) {
    varianceRZ = sigmaY * sigmaY;
  } else if (type == GbtsLayerType::Endcap) {
    varianceRZ = sigmaY * ts.y[1];
    varianceRZ = varianceRZ * varianceRZ;
  } else {
    throw std::runtime_error("invalid layer type");
  }

  const float invResidualVarianceX =
      1.0f / (predictedCovX[0][0] + sigmaX * sigmaX);

  const float invResidualVarianceY = 1.0f / (predictedCovY[0][0] + varianceRZ);

  const float dChi2X = residualX * residualX * invResidualVarianceX;
  const float dChi2Y = residualY * residualY * invResidualVarianceY;

  if (dChi2X > m_cfg.maxDChi2X || dChi2Y > m_cfg.maxDChi2Y) {
    return false;
  }

  // The path search in GraphBasedTrackSeeder::extractBestPathPerEdge skips
  // paths that cannot beat the best one found, assuming a hit adds at most
  // addHit here. A score that can add more must change or remove that bound.
  ts.j += m_cfg.addHit - dChi2X * m_cfg.weightX - dChi2Y * m_cfg.weightY;

  // state update

  const std::array<float, 3> gainX = {
      invResidualVarianceX * predictedCovX[0][0],
      invResidualVarianceX * predictedCovX[0][1],
      invResidualVarianceX * predictedCovX[0][2]};
  const std::array<float, 2> gainY = {
      invResidualVarianceY * predictedCovY[0][0],
      invResidualVarianceY * predictedCovY[0][1]};

  for (std::uint32_t i = 0; i < 3; ++i) {
    ts.x[i] = predictedX[i] + gainX[i] * residualX;
  }

  if (std::abs(ts.x[2]) > m_cfg.maxCurvature) {
    return false;
  }

  for (std::uint32_t i = 0; i < 2; ++i) {
    ts.y[i] = predictedY[i] + gainY[i] * residualY;
  }

  const float z0 = ts.y[0] - refY * ts.y[1];

  const GbtsZ0Range& z0Range = m_geometry->z0Range();
  if (z0 < z0Range.min || z0 > z0Range.max) {
    return false;
  }

  for (std::uint32_t i = 0; i < 3; ++i) {
    for (std::uint32_t j = 0; j < 3; ++j) {
      ts.cx[i][j] = predictedCovX[i][j] - gainX[i] * measurementCovX[j];
    }
  }

  for (std::uint32_t i = 0; i < 2; ++i) {
    for (std::uint32_t j = 0; j < 2; ++j) {
      ts.cy[i][j] = predictedCovY[i][j] - gainY[i] * measurementCovY[j];
    }
  }
  ts.refX = refX;
  ts.refY = refY;
  return true;
}

GbtsLayerType GbtsTrackingFilter::getLayerType(
    const GbtsLayerIndex layerIndex) const {
  return m_geometry->layerDescription(layerIndex).type;
}

}  // namespace Acts::Experimental
