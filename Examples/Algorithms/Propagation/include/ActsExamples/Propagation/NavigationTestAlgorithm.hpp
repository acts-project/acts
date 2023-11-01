// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ActsExamples {

/// @brief this test algorithm performs test propagation
/// within the Acts::Propagator
///
/// If the propagator is equipped appropriately, it can
/// also be used to test the Extrapolator within the geomtetry
class NavigationTestAlgorithm : public IAlgorithm {
 public:
  struct Config {
    /// how to set it up
    std::shared_ptr<RandomNumbers> randomNumberSvc = nullptr;

    /// Tracking geometry input to determine extent
    /// @note We require the top level volume to be a cylinder
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;

    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField = nullptr;

    /// number of particles
    size_t ntests = 100;
  };

  /// Constructor
  /// @param [in] config is the configuration struct
  /// @param [in] loglevel is the logging level
  NavigationTestAlgorithm(const Config& config, Acts::Logging::Level level);

  /// Framework execute method
  /// @param [in] the algorithm context for event consistency
  /// @return is a process code indicating success or not
  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;  ///< the config class

  double m_rMin;
  double m_rMax;
  double m_halfLengthZ;
};

}  // namespace ActsExamples
