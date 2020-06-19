// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <optional>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/DebugOutputActor.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts::UnitLiterals;

namespace FW {

/// Using some short hands for Recorded Material
using RecordedMaterial = Acts::MaterialInteractor::result_type;

/// And recorded material track
/// - this is start:  position, start momentum
///   and the Recorded material
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3D, Acts::Vector3D>, RecordedMaterial>;

/// Finally the output of the propagation test
using PropagationOutput =
    std::pair<std::vector<Acts::detail::Step>, RecordedMaterial>;

/// @brief this test algorithm performs test propagation
/// within the Acts::Propagator
///
/// If the propagator is equipped appropriately, it can
/// also be used to test the Extrapolator within the geomtetry
///
/// @tparam propagator_t Type of the Propagator to be tested
template <typename propagator_t>
class PropagationAlgorithm : public BareAlgorithm {
 public:
  struct Config {
    // create a config object with the propagator
    Config(propagator_t prop) : propagator(std::move(prop)) {}

    /// the propagors to be tested
    propagator_t propagator;

    /// how to set it up
    std::shared_ptr<RandomNumbers> randomNumberSvc = nullptr;

    /// proapgation mode
    int mode = 0;
    /// debug output
    bool debugOutput = false;
    /// Modify the behavior of the material interaction: energy loss
    bool energyLoss = false;
    /// Modify the behavior of the material interaction: scattering
    bool multipleScattering = false;
    /// Modify the behavior of the material interaction: record
    bool recordMaterialInteractions = false;

    /// number of particles
    size_t ntests = 100;
    /// d0 gaussian sigma
    double d0Sigma = 15_um;
    /// z0 gaussian sigma
    double z0Sigma = 55_mm;
    /// phi gaussian sigma (used for covariance transport)
    double phiSigma = 0.0001;
    /// theta gaussian sigma (used for covariance transport)
    double thetaSigma = 0.0001;
    /// qp gaussian sigma (used for covariance transport)
    double qpSigma = 0.00001 / 1_GeV;
    /// t gaussian sigma (used for covariance transport)
    double tSigma = 1_ns;
    /// phi range
    std::pair<double, double> phiRange = {-M_PI, M_PI};
    /// eta range
    std::pair<double, double> etaRange = {-4., 4.};
    /// pt range
    std::pair<double, double> ptRange = {100_MeV, 100_GeV};
    /// looper protection
    double ptLoopers = 300_MeV;

    /// Max step size steering
    double maxStepSize = 3_m;

    /// The step collection to be stored
    std::string propagationStepCollection = "PropagationSteps";

    /// The material collection to be stored
    std::string propagationMaterialCollection = "RecordedMaterialTracks";

    /// covariance transport
    bool covarianceTransport = false;

    /// The covariance values
    Acts::BoundVector covariances = Acts::BoundVector::Zero();

    /// The correlation terms
    Acts::BoundSymMatrix correlations = Acts::BoundSymMatrix::Identity();
  };

  /// Constructor
  /// @param [in] cnf is the configuration struct
  /// @param [in] loglevel is the loggin level
  PropagationAlgorithm(const Config& cnf, Acts::Logging::Level loglevel);

  /// Framework execute method
  /// @param [in] the algorithm context for event consistency
  /// @return is a process code indicating succes or not
  FW::ProcessCode execute(const AlgorithmContext& context) const final override;

 private:
  Config m_cfg;  ///< the config class

  /// Private helper method to create a corrleated covariance matrix
  /// @param[in] rnd is the random engine
  /// @param[in] gauss is a gaussian distribution to draw from
  std::optional<Acts::BoundSymMatrix> generateCovariance(
      FW::RandomEngine& rnd, std::normal_distribution<double>& gauss) const;

  /// Templated execute test method for
  /// charged and netural particles
  ///
  // @tparam parameters_t type of the parameters objects (charged/neutra;)
  ///
  /// @param [in] context The Context for this call
  /// @param [in] startParameters the start parameters
  /// @param [in] pathLengthe the path limit of this propagation
  ///
  /// @return collection of Propagation steps for further analysis
  template <typename parameters_t>
  PropagationOutput executeTest(
      const AlgorithmContext& context, const parameters_t& startParameters,
      double pathLength = std::numeric_limits<double>::max()) const;
};

#include "PropagationAlgorithm.ipp"

}  // namespace FW
