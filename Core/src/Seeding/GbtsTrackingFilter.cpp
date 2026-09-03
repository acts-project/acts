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

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <utility>

namespace Acts::Experimental::detail {

void detail::GbtsEdgeState::initialize(const detail::GbtsEdge& pS,
                                       const detail::GbtsNodeView& nodeView) {
  initialized = true;

  j = 0;
  vs.clear();

  // n2->n1

  const detail::GbtsNodeProxy n1 = nodeView[pS.n1];
  const detail::GbtsNodeProxy n2 = nodeView[pS.n2];

  const float dx = n1.x() - n2.x();
  const float dy = n1.y() - n2.y();
  const float L = std::sqrt(dx * dx + dy * dy);

  s = dy / L;
  c = dx / L;

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
  cx[0][0] = 0.25f;
  cx[1][1] = 0.001f;
  cx[2][2] = 0.001f;

  cy = {};
  cy[0][0] = 1.5f;
  cy[1][1] = 0.001f;
}

}  // namespace Acts::Experimental::detail

namespace Acts::Experimental {

GbtsTrackingFilter::GbtsTrackingFilter(
    const Config& config, const std::shared_ptr<const GbtsGeometry>& geometry,
    std::unique_ptr<const Logger> logger)
    : m_cfg(config), m_geometry(geometry), m_logger(std::move(logger)) {}

detail::GbtsEdgeState GbtsTrackingFilter::followTrack(
    State& state, const detail::GbtsNodeView& nodeView,
    std::vector<detail::GbtsEdge>& sb, detail::GbtsEdge& pS) const {
  if (pS.level == -1) {
    // already collected
    return detail::GbtsEdgeState(false);
  }

  state.globalStateCounter = 0;

  // create track state

  detail::GbtsEdgeState& pInitState =
      state.stateStore[state.globalStateCounter];
  ++state.globalStateCounter;

  pInitState.initialize(pS, nodeView);

  state.stateVec.clear();

  // recursive branching and propagation

  propagate(state, nodeView, sb, pS, pInitState);

  if (state.stateVec.empty()) {
    return detail::GbtsEdgeState(false);
  }

  std::ranges::sort(state.stateVec, std::ranges::greater{},
                    [](const detail::GbtsEdgeState* s) { return s->j; });

  state.globalStateCounter = 0;

  return *state.stateVec.front();
}

void GbtsTrackingFilter::propagate(State& state,
                                   const detail::GbtsNodeView& nodeView,
                                   std::vector<detail::GbtsEdge>& sb,
                                   detail::GbtsEdge& pS,
                                   detail::GbtsEdgeState& ts) const {
  if (state.globalStateCounter >= detail::kGbtsMaxEdgeStates) {
    return;
  }

  detail::GbtsEdgeState& newTs = state.stateStore[state.globalStateCounter];
  ++state.globalStateCounter;
  newTs = ts;

  newTs.vs.push_back(&pS);

  // update using n1 of the segment
  bool accepted = update(nodeView, pS, newTs);

  if (!accepted) {
    // stop further propagation
    return;
  }

  const std::int32_t level = pS.level;

  std::vector<detail::GbtsEdge*> lCont;

  // loop over the neighbours of this segment
  for (std::uint32_t nIdx = 0; nIdx < pS.nNei; ++nIdx) {
    const std::uint32_t nextSegmentIdx = pS.vNei[nIdx];

    detail::GbtsEdge& pN = sb[nextSegmentIdx];

    if (pN.level == -1) {
      // already collected
      continue;
    }

    if (pN.level == level - 1) {
      lCont.push_back(&pN);
    }
  }

  // the end of chain
  if (lCont.empty()) {
    // store in the vector
    if (state.globalStateCounter < detail::kGbtsMaxEdgeStates) {
      if (state.stateVec.empty()) {
        // add the first segment state
        detail::GbtsEdgeState* p = &state.stateStore[state.globalStateCounter];
        ++state.globalStateCounter;
        *p = newTs;
        state.stateVec.push_back(p);
      } else {
        // compare with the best and add
        const float bestSoFar = state.stateVec.front()->j;
        if (newTs.j > bestSoFar) {
          detail::GbtsEdgeState* p =
              &state.stateStore[state.globalStateCounter];
          ++state.globalStateCounter;
          *p = newTs;
          state.stateVec.push_back(p);
        }
      }
    }
  } else {
    // branching
    for (detail::GbtsEdge* sIt : lCont) {
      // recursive call
      propagate(state, nodeView, sb, *sIt, newTs);
    }
  }
}

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

  std::array<float, 3> X{};
  std::array<float, 2> Y{};
  std::array<std::array<float, 3>, 3> Cx{};
  std::array<std::array<float, 2>, 2> Cy{};

  const float x = n1.x();
  const float y = n1.y();
  const float z = n1.z();
  const float r = n1.r();

  const float refX = x * ts.c + y * ts.s;
  const float mx = -x * ts.s + y * ts.c;  // measured X[0]
  const float refY = r;
  const float my = z;  // measured Y[0]

  const float A = refX - ts.refX;
  const float B = 0.5 * A * A;
  const float dr = refY - ts.refY;

  X[0] = ts.x[0] + ts.x[1] * A + ts.x[2] * B;
  X[1] = ts.x[1] + ts.x[2] * A;
  X[2] = ts.x[2];

  Cx[0][0] = ts.cx[0][0] + 2 * ts.cx[0][1] * A + 2 * ts.cx[0][2] * B +
             A * A * ts.cx[1][1] + 2 * A * B * ts.cx[1][2] +
             B * B * ts.cx[2][2];
  Cx[0][1] = Cx[1][0] = ts.cx[0][1] + ts.cx[1][1] * A + ts.cx[1][2] * B +
                        ts.cx[0][2] * A + A * A * ts.cx[1][2] +
                        A * B * ts.cx[2][2];
  Cx[0][2] = Cx[2][0] = ts.cx[0][2] + ts.cx[1][2] * A + ts.cx[2][2] * B;

  Cx[1][1] = ts.cx[1][1] + 2 * A * ts.cx[1][2] + A * A * ts.cx[2][2];
  Cx[1][2] = Cx[2][1] = ts.cx[1][2] + ts.cx[2][2] * A;

  Cx[2][2] = ts.cx[2][2];

  Y[0] = ts.y[0] + ts.y[1] * dr;
  Y[1] = ts.y[1];

  Cy[0][0] = ts.cy[0][0] + 2 * ts.cy[0][1] * dr + dr * dr * ts.cy[1][1];
  Cy[0][1] = Cy[1][0] = ts.cy[0][1] + dr * ts.cy[1][1];
  Cy[1][1] = ts.cy[1][1];

  // chi2 test

  const float resid_x = mx - X[0];
  const float resid_y = my - Y[0];

  const std::array<float, 3> CHx = {Cx[0][0], Cx[0][1], Cx[0][2]};
  const std::array<float, 2> CHy = {Cy[0][0], Cy[0][1]};

  float sigma_rz = 0;

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
    sigma_rz = sigmaY * sigmaY;
  } else if (type == GbtsLayerType::Endcap) {
    sigma_rz = sigmaY * ts.y[1];
    sigma_rz = sigma_rz * sigma_rz;
  } else {
    throw std::runtime_error("invalid layer type");
  }

  const float Dx = 1.0f / (Cx[0][0] + sigmaX * sigmaX);

  const float Dy = 1.0f / (Cy[0][0] + sigma_rz);

  const float dchi2_x = resid_x * resid_x * Dx;
  const float dchi2_y = resid_y * resid_y * Dy;

  if (dchi2_x > m_cfg.maxDChi2X || dchi2_y > m_cfg.maxDChi2Y) {
    return false;
  }

  ts.j += m_cfg.addHit - dchi2_x * m_cfg.weightX - dchi2_y * m_cfg.weightY;

  // state update

  const std::array<float, 3> Kx = {Dx * Cx[0][0], Dx * Cx[0][1], Dx * Cx[0][2]};
  const std::array<float, 2> Ky = {Dy * Cy[0][0], Dy * Cy[0][1]};

  for (std::uint32_t i = 0; i < 3; ++i) {
    ts.x[i] = X[i] + Kx[i] * resid_x;
  }

  if (std::abs(ts.x[2]) > m_cfg.maxCurvature) {
    return false;
  }

  for (std::uint32_t i = 0; i < 2; ++i) {
    ts.y[i] = Y[i] + Ky[i] * resid_y;
  }

  const float z0 = ts.y[0] - refY * ts.y[1];

  if (std::abs(z0) > m_cfg.maxZ0) {
    return false;
  }

  for (std::uint32_t i = 0; i < 3; ++i) {
    for (std::uint32_t j = 0; j < 3; ++j) {
      ts.cx[i][j] = Cx[i][j] - Kx[i] * CHx[j];
    }
  }

  for (std::uint32_t i = 0; i < 2; ++i) {
    for (std::uint32_t j = 0; j < 2; ++j) {
      ts.cy[i][j] = Cy[i][j] - Ky[i] * CHy[j];
    }
  }
  ts.refX = refX;
  ts.refY = refY;
  return true;
}

GbtsLayerType GbtsTrackingFilter::getLayerType(
    const std::uint32_t layerIndex) const {
  return m_geometry->layerDescriptionByIndex(layerIndex).type;
}

}  // namespace Acts::Experimental
