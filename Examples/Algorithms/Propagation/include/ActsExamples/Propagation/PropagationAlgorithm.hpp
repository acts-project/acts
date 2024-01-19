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

class PropagatorInterface;
struct AlgorithmContext;

/// Using some short hands for Recorded Material
using RecordedMaterial = Acts::MaterialInteractor::result_type;

/// And recorded material track
/// - this is start:  position, start momentum
///   and the Recorded material
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3, Acts::Vector3>, RecordedMaterial>;

/// Finally the output of the propagation test
using PropagationOutput =
    std::pair<std::vector<Acts::detail::Step>, RecordedMaterial>;

/// @brief this test algorithm performs test propagation
/// within the Acts::Propagator
///
/// If the propagator is equipped appropriately, it can
/// also be used to test the Extrapolator within the geomtetry
class PropagationAlgorithm : public IAlgorithm {
 public:
  struct Config {
    /// Instance of a propagator wrapper that performs the actual propagation
    std::shared_ptr<PropagatorInterface> propagatorImpl = nullptr;

    /// how to set it up
    std::shared_ptr<RandomNumbers> randomNumberSvc = nullptr;
    /// proapgation mode
    int mode = 0;
    /// Switch the logger to sterile
    bool sterileLogger = false;
    /// debug output
    bool debugOutput = false;
    /// Modify the behavior of the material interaction: energy loss
    bool energyLoss = true;
    /// Modify the behavior of the material interaction: scattering
    bool multipleScattering = true;
    /// Modify the behavior of the material interaction: record
    bool recordMaterialInteractions = true;

    /// number of particles
    std::size_t ntests = 100;
    /// d0 gaussian sigma
    double d0Sigma = 15 * Acts::UnitConstants::um;
    /// z0 gaussian sigma
    double z0Sigma = 55 * Acts::UnitConstants::mm;
    /// phi gaussian sigma (used for covariance transport)
    double phiSigma = 0.001;
    /// theta gaussian sigma (used for covariance transport)
    double thetaSigma = 0.001;
    /// qp gaussian sigma (used for covariance transport)
    double qpSigma = 0.0001 / 1 * Acts::UnitConstants::GeV;
    /// t gaussian sigma (used for covariance transport)
    double tSigma = 1 * Acts::UnitConstants::ns;
    /// phi range
    std::pair<double, double> phiRange = {-M_PI, M_PI};
    /// eta range
    std::pair<double, double> etaRange = {-4., 4.};
    /// pt range
    std::pair<double, double> ptRange = {100 * Acts::UnitConstants::MeV,
                                         100 * Acts::UnitConstants::GeV};
    /// particle hypothesis
    Acts::ParticleHypothesis particleHypothesis =
        Acts::ParticleHypothesis::pion();
    /// looper protection
    double ptLoopers = 500 * Acts::UnitConstants::MeV;

    /// Max step size steering
    double maxStepSize = 3 * Acts::UnitConstants::m;

    /// The step collection to be stored
    std::string propagationStepCollection = "PropagationSteps";

    /// The material collection to be stored
    std::string propagationMaterialCollection = "RecordedMaterialTracks";

    /// covariance transport
    bool covarianceTransport = false;

    /// The covariance values
    Acts::BoundVector covariances = Acts::BoundVector::Zero();

    /// The correlation terms
    Acts::BoundSquareMatrix correlations = Acts::BoundSquareMatrix::Identity();
  };

  /// Constructor
  /// @param [in] config is the configuration struct
  /// @param [in] loglevel is the logging level
  PropagationAlgorithm(const Config& config, Acts::Logging::Level level);

  /// Framework execute method
  /// @param [in] the algorithm context for event consistency
  /// @return is a process code indicating success or not
  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;  ///< the config class

  WriteDataHandle<std::vector<std::vector<Acts::detail::Step>>>
      m_outpoutPropagationSteps{this, "OutputPropagationSteps"};

  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_recordedMaterial{this, "RecordedMaterial"};

  /// Private helper method to create a corrleated covariance matrix
  /// @param[in] rnd is the random engine
  /// @param[in] gauss is a gaussian distribution to draw from
  std::optional<Acts::BoundSquareMatrix> generateCovariance(
      ActsExamples::RandomEngine& rnd,
      std::normal_distribution<double>& gauss) const;
};

}  // namespace ActsExamples
