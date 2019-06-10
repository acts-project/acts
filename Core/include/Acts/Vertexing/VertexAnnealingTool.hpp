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

/// @brief Implements a thermodynamic annealing scheme for the track
///   weight factors used in the MultiAdaptiveVertexFitter in such a way
///   that with high temperature values (at the beginning) only a slight
///   preference is given to tracks compatible with the estimated vertex
///   position. With lower temperatures the weighting get stricter such
///   that all incompatible tracks will be dropped at the end while
///   keeping all compatible tracks with a track weight of 1.
///   Ref. (1): CERN-THESIS-2010-027, Author: Piacquadio, Giacinto:
///   `Identification of b-jets and investigation of the discovery potential
///   of a Higgs boson in the WH−−>lvbb¯ channel with the ATLAS experiment`
class VertexAnnealingTool {
 public:
  struct State {
    /// Points to current temperature value in m_cfg.setOfTemperatures
    unsigned int currentTemperatureIndex = 0;

    /// Checks if equilibrium is reached
    bool equilibriumReached = false;
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
  /// @param chi2 Chi^2 for current track, i.e. compatibility
  /// of track to current vertex candidate
  /// @param allChi2 Vector of all chi^2 values, i.e. compatibilities
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

  /// @brief Gaussian function for weight calculation
  ///
  /// @param chi2 Chi^2 value
  /// @param temp Temperature value
  ///
  /// @return exp(-1./2. * chi2 / temp)
  double gaussFunc(const double chi2, const double temp) const {
    return std::exp(-1. / 2. * chi2 / temp);
  }
};
