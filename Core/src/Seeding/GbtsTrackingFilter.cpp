// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsTrackingFilter.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <list>

namespace Acts::Experimental {

void GbtsEdgeState::initialize(GbtsEdge& pS) {
  m_initialized = true;

  m_J = 0.0;
  m_vs.clear();

  // n2->n1

  float dx = pS.m_n1->x() - pS.m_n2->x();
  float dy = pS.m_n1->y() - pS.m_n2->y();
  float L = std::sqrt(dx * dx + dy * dy);

  m_s = dy / L;
  m_c = dx / L;

  // transform for extrapolation and update
  //  x' =  x*m_c + y*m_s
  //  y' = -x*m_s + y*m_c

  m_refY = pS.m_n2->r();
  m_refX = pS.m_n2->x() * m_c + pS.m_n2->y() * m_s;

  // X-state: y, dy/dx, d2y/dx2

  m_X[0] = -pS.m_n2->x() * m_s + pS.m_n2->y() * m_c;
  m_X[1] = 0.0;
  m_X[2] = 0.0;

  // Y-state: z, dz/dr

  m_Y[0] = pS.m_n2->z();
  m_Y[1] = (pS.m_n1->z() - pS.m_n2->z()) / (pS.m_n1->r() - pS.m_n2->r());

  memset(&m_Cx[0][0], 0, sizeof(m_Cx));
  memset(&m_Cy[0][0], 0, sizeof(m_Cy));

  m_Cx[0][0] = 0.25;
  m_Cx[1][1] = 0.001;
  m_Cx[2][2] = 0.001;

  m_Cy[0][0] = 1.5;
  m_Cy[1][1] = 0.001;
}

void GbtsEdgeState::clone(const GbtsEdgeState& st) {
  memcpy(&m_X[0], &st.m_X[0], sizeof(m_X));
  memcpy(&m_Y[0], &st.m_Y[0], sizeof(m_Y));
  memcpy(&m_Cx[0][0], &st.m_Cx[0][0], sizeof(m_Cx));
  memcpy(&m_Cy[0][0], &st.m_Cy[0][0], sizeof(m_Cy));
  m_refX = st.m_refX;
  m_refY = st.m_refY;
  m_c = st.m_c;
  m_s = st.m_s;
  m_J = st.m_J;
  m_vs.clear();
  m_vs.reserve(st.m_vs.size());
  std::copy(st.m_vs.begin(), st.m_vs.end(), std::back_inserter(m_vs));

  m_initialized = true;
}

GbtsTrackingFilter::GbtsTrackingFilter(const std::vector<TrigInDetSiLayer>& g,
                                       std::vector<GbtsEdge>& sb,
                                       const SeedFinderGbtsConfig& config)
    : m_geo(g), m_segStore(sb), m_config(config) {}

void GbtsTrackingFilter::followTrack(GbtsEdge& pS, GbtsEdgeState& output) {
  if (pS.m_level == -1) {
    return;  // already collected
  }

  m_globalStateCounter = 0;

  // create track state

  GbtsEdgeState& pInitState = m_stateStore[m_globalStateCounter++];

  pInitState.initialize(pS);

  m_stateVec.clear();

  // recursive branching and propagation

  propagate(pS, pInitState);

  if (m_stateVec.empty()) {
    return;
  }

  std::sort(m_stateVec.begin(), m_stateVec.end(),
            typename GbtsEdgeState::Compare());

  GbtsEdgeState* best = (*m_stateVec.begin());

  output.clone(*best);

  m_globalStateCounter = 0;
}

void GbtsTrackingFilter::propagate(GbtsEdge& pS, GbtsEdgeState& ts) {
  if (m_globalStateCounter >= GbtsMaxEdgeState) {
    return;
  }

  GbtsEdgeState* p_new_ts = &m_stateStore[m_globalStateCounter++];

  GbtsEdgeState& new_ts = *p_new_ts;
  new_ts.clone(ts);

  new_ts.m_vs.push_back(&pS);

  bool accepted = update(pS, new_ts);  // update using n1 of the segment

  if (!accepted) {
    return;  // stop further propagation
  }

  std::int32_t level = pS.m_level;

  std::list<GbtsEdge*> lCont;

  // loop over the neighbours of this segment
  for (std::uint32_t nIdx = 0; nIdx < pS.m_nNei; ++nIdx) {
    std::uint32_t nextSegmentIdx = pS.m_vNei[nIdx];

    GbtsEdge* pN = &(m_segStore[nextSegmentIdx]);

    if (pN->m_level == -1) {
      continue;  // already collected
    }

    if (pN->m_level == level - 1) {
      lCont.push_back(pN);
    }
  }

  if (lCont.empty()) {  // the end of chain

    // store in the vector
    if (m_globalStateCounter < GbtsMaxEdgeState) {
      if (m_stateVec.empty()) {  // add the first segment state
        GbtsEdgeState* p = &m_stateStore[m_globalStateCounter++];
        p->clone(new_ts);
        m_stateVec.push_back(p);
      } else {  // compare with the best and add
        float best_so_far = (*m_stateVec.begin())->m_J;
        if (new_ts.m_J > best_so_far) {
          GbtsEdgeState* p = &m_stateStore[m_globalStateCounter++];
          p->clone(new_ts);
          m_stateVec.push_back(p);
        }
      }
    }
  } else {  // branching

    for (const auto sIt : lCont) {
      propagate(*sIt, new_ts);  // recursive call
    }
  }
}

bool GbtsTrackingFilter::update(GbtsEdge& pS, GbtsEdgeState& ts) {
  if (ts.m_Cx[2][2] < 0.0 || ts.m_Cx[1][1] < 0.0 || ts.m_Cx[0][0] < 0.0) {
    std::cout << "Negative cov_x" << std::endl;
  }

  if (ts.m_Cy[1][1] < 0.0 || ts.m_Cy[0][0] < 0.0) {
    std::cout << "Negative cov_y" << std::endl;
  }

  // add ms.

  float tau2 = ts.m_Y[1] * ts.m_Y[1];
  float invSin2 = 1 + tau2;

  std::int32_t type1 = getLayerType(pS.m_n2->layer());  // 0 - barrel

  float lenCorr = type1 == 0 ? invSin2 : invSin2 / tau2;

  float minPtFrac = std::abs(ts.m_X[2]) / m_config.max_curvature;

  float corrMS = m_config.sigmaMS * minPtFrac;

  float sigma2 = m_config.radLen * lenCorr * corrMS * corrMS;  // /invSin2;

  ts.m_Cx[1][1] += sigma2;

  ts.m_Cy[1][1] += sigma2;

  // extrapolation

  float X[3], Y[2];
  float Cx[3][3], Cy[2][2];

  float refX{}, refY{}, mx{}, my{};

  float x{}, y{}, z{}, r{};

  x = pS.m_n1->x();
  y = pS.m_n1->y();
  z = pS.m_n1->z();
  r = pS.m_n1->r();

  refX = x * ts.m_c + y * ts.m_s;
  mx = -x * ts.m_s + y * ts.m_c;  // measured X[0]
  refY = r;
  my = z;  // measured Y[0]

  float A = refX - ts.m_refX;
  float B = 0.5 * A * A;
  float dr = refY - ts.m_refY;

  X[0] = ts.m_X[0] + ts.m_X[1] * A + ts.m_X[2] * B;
  X[1] = ts.m_X[1] + ts.m_X[2] * A;
  X[2] = ts.m_X[2];

  Cx[0][0] = ts.m_Cx[0][0] + 2 * ts.m_Cx[0][1] * A + 2 * ts.m_Cx[0][2] * B +
             A * A * ts.m_Cx[1][1] + 2 * A * B * ts.m_Cx[1][2] +
             B * B * ts.m_Cx[2][2];
  Cx[0][1] = Cx[1][0] = ts.m_Cx[0][1] + ts.m_Cx[1][1] * A + ts.m_Cx[1][2] * B +
                        ts.m_Cx[0][2] * A + A * A * ts.m_Cx[1][2] +
                        A * B * ts.m_Cx[2][2];
  Cx[0][2] = Cx[2][0] = ts.m_Cx[0][2] + ts.m_Cx[1][2] * A + ts.m_Cx[2][2] * B;

  Cx[1][1] = ts.m_Cx[1][1] + 2 * A * ts.m_Cx[1][2] + A * A * ts.m_Cx[2][2];
  Cx[1][2] = Cx[2][1] = ts.m_Cx[1][2] + ts.m_Cx[2][2] * A;

  Cx[2][2] = ts.m_Cx[2][2];

  Y[0] = ts.m_Y[0] + ts.m_Y[1] * dr;
  Y[1] = ts.m_Y[1];

  Cy[0][0] = ts.m_Cy[0][0] + 2 * ts.m_Cy[0][1] * dr + dr * dr * ts.m_Cy[1][1];
  Cy[0][1] = Cy[1][0] = ts.m_Cy[0][1] + dr * ts.m_Cy[1][1];
  Cy[1][1] = ts.m_Cy[1][1];

  // chi2 test

  float resid_x = mx - X[0];
  float resid_y = my - Y[0];

  float CHx[3] = {Cx[0][0], Cx[0][1], Cx[0][2]};
  float CHy[2] = {Cy[0][0], Cy[0][1]};

  float sigma_rz = 0.0;

  std::int32_t type = getLayerType(pS.m_n1->layer());

  if (type == 0) {  // barrel TO-DO: split into barrel Pixel and barrel SCT
    sigma_rz = m_config.sigma_y * m_config.sigma_y;
  } else {
    sigma_rz = m_config.sigma_y * ts.m_Y[1];
    sigma_rz = sigma_rz * sigma_rz;
  }

  float Dx = 1.0 / (Cx[0][0] + m_config.sigma_x * m_config.sigma_x);

  float Dy = 1.0 / (Cy[0][0] + sigma_rz);

  float dchi2_x = resid_x * resid_x * Dx;
  float dchi2_y = resid_y * resid_y * Dy;

  if (dchi2_x > m_config.maxDChi2_x || dchi2_y > m_config.maxDChi2_y) {
    return false;
  }

  ts.m_J += m_config.add_hit - dchi2_x * m_config.weight_x -
            dchi2_y * m_config.weight_y;

  // state update

  float Kx[3] = {Dx * Cx[0][0], Dx * Cx[0][1], Dx * Cx[0][2]};
  float Ky[2] = {Dy * Cy[0][0], Dy * Cy[0][1]};

  for (std::uint32_t i = 0; i < 3; ++i) {
    ts.m_X[i] = X[i] + Kx[i] * resid_x;
  }

  if (std::abs(ts.m_X[2]) > m_config.max_curvature) {
    return false;
  }

  for (std::uint32_t i = 0; i < 2; ++i) {
    ts.m_Y[i] = Y[i] + Ky[i] * resid_y;
  }

  float z0 = ts.m_Y[0] - refY * ts.m_Y[1];

  if (std::abs(z0) > m_config.max_z0) {
    return false;
  }

  for (std::uint32_t i = 0; i < 3; ++i) {
    for (std::uint32_t j = 0; j < 3; ++j) {
      ts.m_Cx[i][j] = Cx[i][j] - Kx[i] * CHx[j];
    }
  }

  for (std::uint32_t i = 0; i < 2; ++i) {
    for (std::uint32_t j = 0; j < 2; ++j) {
      ts.m_Cy[i][j] = Cy[i][j] - Ky[i] * CHy[j];
    }
  }
  ts.m_refX = refX;
  ts.m_refY = refY;
  return true;
}

std::uint32_t GbtsTrackingFilter::getLayerType(std::uint32_t l) {
  return m_geo.at(l).m_type;
}

}  // namespace Acts::Experimental
