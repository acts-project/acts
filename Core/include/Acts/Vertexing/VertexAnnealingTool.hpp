// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include "Acts/Utilities/Definitions.hpp"

// TODO remove
#include <iostream>

class VertexAnnealingTool {
 public:
  struct State {
    /// Constructor
    State() : currentTemperatureIndex(0), equilibriumReached(false) {}

    /// Points to current temperature value in m_cfg.setOfTemperatures
    unsigned int currentTemperatureIndex;

    /// Checks if equilibrium is reached
    bool equilibriumReached;
  };

  struct Config {
    // Config constructor with default temperature list: {64.,16.,4.,2.,1.5,1.}
    Config(const std::vector<double>& temperatures = {64., 16., 4., 2., 1.5,
                                                      1.})
        : setOfTemperatures(temperatures) {
      if (setOfTemperatures.empty()) {
        // TODO: error handling, error msg
        std::cout << "ERROR: setOfTemperatures must not be empty" << std::endl;
      }
    }

    /// Insensitivity of calculated weight at cutoff
    double cutOff = 9.;

    /// Set of temperatures, annealing starts at setOfTemperatures[0]
    /// and anneals towards setOfTemperatures[last]
    std::vector<double> setOfTemperatures;
  };

  /// Constructor
  VertexAnnealingTool(const Config& cfg = Config()) : m_cfg(cfg) {}

  /// Resets the annealing process
  void reset(State& state) const;

  /// Does the actual annealing step
  ///
  void anneal(State& state) const;

  /// @brief Weight access
  ///
  /// @param chi2 Chi^2
  /// @param allChi2 Vector of all chi^2
  ///
  /// @return Calculated weight
  double getWeight(State& state, double chi2,
                   const std::vector<double>& allChi2) const;

  /// @brief Weight access
  ///
  /// @param chi2 Chi^2
  ///
  /// @return Calculated weight
  double getWeight(State& state, double chi2) const;

 private:
  /// Configuration object
  Config m_cfg;

  /// @brief Gaussian function for weight calculation
  ///
  /// @param val1 Value 1
  /// @param val2 Value 2
  ///
  /// @return exp(-1./2. * val1 / val2)
  double gaussFunc(const double val1, const double val2) const {
    return std::exp(-1. / 2. * val1 / val2);
  }
};
