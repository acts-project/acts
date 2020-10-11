// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/AnnealingUtility.hpp"

/// @brief Gaussian-like function for weight calculation
/// Note: Factor 2 in denominator is included in temperature
///
/// @param chi2 Chi^2 value
/// @param temp Temperature value
///
/// @return exp(-chi2 / temp)
static double gaussFunc(double chi2, double temp) {
  return std::exp(-chi2 / temp);
}

void Acts::AnnealingUtility::anneal(State& state) const {
  if (state.currentTemperatureIndex < m_cfg.setOfTemperatures.size() - 1) {
    ++state.currentTemperatureIndex;
  } else {
    state.equilibriumReached = true;
  }
}

double Acts::AnnealingUtility::getWeight(
    State& state, double chi2, const std::vector<double>& allChi2) const {
  unsigned int idx = state.currentTemperatureIndex;
  const double currentTemp = 2. * m_cfg.setOfTemperatures[idx];

  double num = gaussFunc(chi2, currentTemp);

  double denom = m_cfg.gaussCutTempVec[idx];

  for (double val : allChi2) {
    denom += gaussFunc(val, currentTemp);
  }

  return num / denom;
}

double Acts::AnnealingUtility::getWeight(State& state, double chi2) const {
  const double currentTemp =
      m_cfg.setOfTemperatures[state.currentTemperatureIndex];

  return 1. / (1. + gaussFunc(m_cfg.cutOff - chi2, 2. * currentTemp));
}
