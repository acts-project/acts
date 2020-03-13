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

/// @brief Implements a deterministic thermodynamic annealing scheme
/// Ref. (1): CERN-THESIS-2010-027
class AnnealingUtility {
 public:
  /// @brief The annealing state
  /// Resetting the state is done by just creating a new instance
  struct State {
    // Points to current temperature value in m_cfg.setOfTemperatures
    unsigned int currentTemperatureIndex{0};

    // Checks if equilibrium is reached
    bool equilibriumReached{false};
  };

  /// @brief The configuration struct
  struct Config {
    // Config constructor with default temperature list: {64.,16.,4.,2.,1.5,1.}
    Config(const std::vector<double>& temperatures = {64., 16., 4., 2., 1.5,
                                                      1.})
        : setOfTemperatures(temperatures) {}

    // Insensitivity of calculated weight at cutoff
    double cutOff{9.};

    // Set of temperatures, annealing starts at setOfTemperatures[0]
    // and anneals towards setOfTemperatures[last]
    std::vector<double> setOfTemperatures;
  };

  /// Constructor
  AnnealingUtility(const Config& cfg = Config()) : m_cfg(cfg) {}

  /// Does the actual annealing step
  void anneal(State& state) const;

  /// @brief Weight access
  ///
  /// @param chi2 Chi^2 for e.g. current track, i.e. compatibility
  /// of track to current vertex candidate
  /// @param allChi2 Vector of all chi^2 values, i.e. e.g. compatibilities
  /// of current track to all vertices it is currently attached to
  ///
  /// @return Calculated weight according to Eq.(5.46) in Ref.(1)
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
};
