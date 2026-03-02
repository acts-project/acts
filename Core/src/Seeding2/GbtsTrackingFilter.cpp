// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/GbtsTrackingFilter.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

namespace Acts::Experimental {

void GbtsEdgeState::initialize(const GbtsEdge& pS) {
  initialized = true;

  j = 0;
  vs.clear();

  // n2->n1

  const float dx = pS.n1->x - pS.n2->x;
  const float dy = pS.n1->y - pS.n2->y;
  const float L = std::sqrt(dx * dx + dy * dy);

  s = dy / L;
  c = dx / L;

  // transform for extrapolation and update
  //  x' =  x*c + y*s
  //  y' = -x*s + y*c

  refY = pS.n2->r;
  refX = pS.n2->x * c + pS.n2->y * s;

  // X-state: y, dy/dx, d2y/dx2

  x[0] = -pS.n2->x * s + pS.n2->y * c;
  x[1] = 0;
  x[2] = 0;

  // Y-state: z, dz/dr

  y[0] = pS.n2->z;
  y[1] = (pS.n1->z - pS.n2->z) / (pS.n1->r - pS.n2->r);

  cx = {};
  cx[0][0] = 0.25f;
  cx[1][1] = 0.001f;
  cx[2][2] = 0.001f;

  cy = {};
  cy[0][0] = 1.5f;
  cy[1][1] = 0.001f;
}

GbtsTrackingFilter::GbtsTrackingFilter(
    const GbtsConfig& config,
    const std::vector<TrigInDetSiLayer>& layerGeometry,
    std::vector<GbtsEdge>& sb)
    : m_cfg(&config), m_layerGeometry(&layerGeometry), m_segStore(&sb) {}

GbtsEdgeState GbtsTrackingFilter::followTrack(GbtsEdge& pS) {
  if (pS.level == -1) {
    // already collected
    return GbtsEdgeState(false);
  }

  m_globalStateCounter = 0;

  // create track state

  GbtsEdgeState& pInitState = m_stateStore[m_globalStateCounter];
  ++m_globalStateCounter;

  pInitState.initialize(pS);

  m_stateVec.clear();

  // recursive branching and propagation

  propagate(pS, pInitState);

  if (m_stateVec.empty()) {
    return GbtsEdgeState(false);
  }

  std::ranges::sort(m_stateVec, std::ranges::greater{},
                    [](const GbtsEdgeState* s) { return s->j; });

  m_globalStateCounter = 0;

  return *m_stateVec.front();
}

void GbtsTrackingFilter::propagate(GbtsEdge& pS, GbtsEdgeState& ts) {
  if (m_globalStateCounter >= GbtsMaxEdgeState) {
    return;
  }

  GbtsEdgeState& newTs = m_stateStore[m_globalStateCounter];
  ++m_globalStateCounter;
  newTs = ts;

  newTs.vs.push_back(&pS);

  // update using n1 of the segment
  bool accepted = update(pS, newTs);

  if (!accepted) {
    // stop further propagation
    return;
  }

  const std::int32_t level = pS.level;

  std::vector<GbtsEdge*> lCont;

  // loop over the neighbours of this segment
  for (std::uint32_t nIdx = 0; nIdx < pS.nNei; ++nIdx) {
    const std::uint32_t nextSegmentIdx = pS.vNei[nIdx];

    GbtsEdge& pN = (*m_segStore)[nextSegmentIdx];

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
    if (m_globalStateCounter < GbtsMaxEdgeState) {
      if (m_stateVec.empty()) {
        // add the first segment state
        GbtsEdgeState* p = &m_stateStore[m_globalStateCounter];
        ++m_globalStateCounter;
        *p = newTs;
        m_stateVec.push_back(p);
      } else {
        // compare with the best and add
        const float bestSoFar = m_stateVec.front()->j;
        if (newTs.j > bestSoFar) {
          GbtsEdgeState* p = &m_stateStore[m_globalStateCounter];
          ++m_globalStateCounter;
          *p = newTs;
          m_stateVec.push_back(p);
        }
      }
    }
  } else {
    // branching
    for (GbtsEdge* sIt : lCont) {
      // recursive call
      propagate(*sIt, newTs);
    }
  }
}

bool GbtsTrackingFilter::update(const GbtsEdge& pS, GbtsEdgeState& ts) const {
  if (ts.cx[2][2] < 0 || ts.cx[1][1] < 0 || ts.cx[0][0] < 0) {
    std::cout << "Negative cov_x" << std::endl;
  }

  if (ts.cy[1][1] < 0 || ts.cy[0][0] < 0) {
    std::cout << "Negative cov_y" << std::endl;
  }

  // add ms.

  const float tau2 = ts.y[1] * ts.y[1];
  const float invSin2 = 1 + tau2;

  const std::int32_t type1 = getLayerType(pS.n2->layer);  // 0 - barrel

  const float lenCorr = type1 == 0 ? invSin2 : invSin2 / tau2;

  const float minPtFrac = std::abs(ts.x[2]) / m_cfg->maxCurvature;

  const float corrMS = m_cfg->sigmaMS * minPtFrac;

  const float sigma2 = m_cfg->radLen * lenCorr * corrMS * corrMS;  // /invSin2

  ts.cx[1][1] += sigma2;

  ts.cy[1][1] += sigma2;

  // extrapolation

  std::array<float, 3> X{};
  std::array<float, 2> Y{};
  std::array<std::array<float, 3>, 3> Cx{};
  std::array<std::array<float, 2>, 2> Cy{};

  const float x = pS.n1->x;
  const float y = pS.n1->y;
  const float z = pS.n1->z;
  const float r = pS.n1->r;

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

  const std::int32_t type = getLayerType(pS.n1->layer);

  if (type == 0) {
    // barrel TODO: split into barrel Pixel and barrel SCT
    sigma_rz = m_cfg->sigmaY * m_cfg->sigmaY;
  } else {
    sigma_rz = m_cfg->sigmaY * ts.y[1];
    sigma_rz = sigma_rz * sigma_rz;
  }

  const float Dx = 1.0 / (Cx[0][0] + m_cfg->sigmaX * m_cfg->sigmaX);

  const float Dy = 1.0 / (Cy[0][0] + sigma_rz);

  const float dchi2_x = resid_x * resid_x * Dx;
  const float dchi2_y = resid_y * resid_y * Dy;

  if (dchi2_x > m_cfg->maxDChi2X || dchi2_y > m_cfg->maxDChi2Y) {
    return false;
  }

  ts.j += m_cfg->addHit - dchi2_x * m_cfg->weightX - dchi2_y * m_cfg->weightY;

  // state update

  const std::array<float, 3> Kx = {Dx * Cx[0][0], Dx * Cx[0][1], Dx * Cx[0][2]};
  const std::array<float, 2> Ky = {Dy * Cy[0][0], Dy * Cy[0][1]};

  for (std::uint32_t i = 0; i < 3; ++i) {
    ts.x[i] = X[i] + Kx[i] * resid_x;
  }

  if (std::abs(ts.x[2]) > m_cfg->maxCurvature) {
    return false;
  }

  for (std::uint32_t i = 0; i < 2; ++i) {
    ts.y[i] = Y[i] + Ky[i] * resid_y;
  }

  const float z0 = ts.y[0] - refY * ts.y[1];

  if (std::abs(z0) > m_cfg->maxZ0) {
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

std::uint32_t GbtsTrackingFilter::getLayerType(const std::uint32_t l) const {
  return m_layerGeometry->at(l).type;
}

}  // namespace Acts::Experimental
