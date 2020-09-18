// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <memory>
#include <string>
#include <unordered_map>

namespace ActsExamples {

/// @brief Digitization Algorithm that turns simulated
/// hits into measuremetns for Fitting
class SmearingAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input collection of simulated hits
    std::string inputSimulatedHits;
    /// Output collection of measuremetns
    std::string outputMeasurements;
    /// Random numbers tool.
    std::shared_ptr<const RandomNumbers> randomNumbers;
  };

  /// Construct the digitization algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SmearingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Build measurement from simulation hits at input
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const final override;

 private:
  Config m_cfg;

  ActsFatras::UncorrelatedHitSmearer m_smearer;

  template <typename... arguments_t>
  auto pixelParSet(arguments_t &&... args) const {
    return m_smearer.smearedParameterSet<Acts::eBoundLoc0, Acts::eBoundLoc1>(
        std::forward<arguments_t>(args)...);
  }

  template <typename... arguments_t>
  auto pixelTimedParSet(arguments_t &&... args) const {
    return m_smearer.smearedParameterSet<Acts::eBoundLoc0, Acts::eBoundLoc1,
                                         Acts::eBoundTime>(
        std::forward<arguments_t>(args)...);
  }

  template <typename... arguments_t>
  auto stripLoc0ParSet(arguments_t &&... args) const {
    return m_smearer.smearedParameterSet<Acts::eBoundLoc0>(
        std::forward<arguments_t>(args)...);
  }

  template <typename... arguments_t>
  auto stripLoc1ParSet(arguments_t &&... args) const {
    return m_smearer.smearedParameterSet<Acts::eBoundLoc1>(
        std::forward<arguments_t>(args)...);
  }

};

}  // namespace ActsExamples