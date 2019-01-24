// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/VertexAnnealingTool.hpp"

void VertexAnnealingTool::reset(State& state) const {
  state.currentTemperatureIndex = 0;
  state.equilibriumReached = false;
}

void VertexAnnealingTool::anneal(State& state) const {
  if (state.currentTemperatureIndex < m_cfg.setOfTemperatures.size() - 1) {
    ++state.currentTemperatureIndex;
  } else {
    state.equilibriumReached = true;
  }
}

double VertexAnnealingTool::getWeight(
    State& state, double chi2, const std::vector<double>& allChi2) const {
  const double currentTemp =
      m_cfg.setOfTemperatures[state.currentTemperatureIndex];

  double allWeights = 0.;
  for (auto val : allChi2) {
    allWeights += gaussFunc(val, currentTemp);
  }

  double actualWeight = gaussFunc(chi2, currentTemp);

  return actualWeight / (gaussFunc(m_cfg.cutOff, currentTemp) + allWeights);
}

double VertexAnnealingTool::getWeight(State& state, double chi2) const {
  const double currentTemp =
      m_cfg.setOfTemperatures[state.currentTemperatureIndex];

  return gaussFunc(chi2, currentTemp) /
         (gaussFunc(m_cfg.cutOff, currentTemp) + gaussFunc(chi2, currentTemp));
}
