// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <iostream>

#include <boost/program_options.hpp>

#include "PropagationAlgorithm.hpp"

namespace ActsExamples {

namespace Options {

/// @brief PropagationAlgorithm options
///
/// @tparam aopt_t Type of the options class from boost
inline void addPropagationOptions(
    boost::program_options::options_description& opt) {
  namespace po = boost::program_options;
  using namespace Acts::UnitLiterals;
  opt.add_options()(
      "prop-debug", po::value<bool>()->default_value(false),
      "Run in debug mode, will create propagation screen output.")(
      "prop-step-collection",
      po::value<std::string>()->default_value("propagation-steps"),
      "Propgation step collection.")(
      "prop-stepper", po::value<int>()->default_value(1),
      "Propgation type: 0 (StraightLine), 1 (Eigen), 2 (Atlas).")(
      "prop-mode", po::value<int>()->default_value(0),
      "Propgation modes: 0 (inside-out), 1 (surface to surface).")(
      "prop-cov", po::value<bool>()->default_value(false),
      "Propagate (random) test covariances.")(
      "prop-energyloss", po::value<bool>()->default_value(true),
      "Apply energy loss correction - in extrapolation mode only.")(
      "prop-scattering", po::value<bool>()->default_value(true),
      "Apply scattering correction - in extrapolation mode only.")(
      "prop-record-material", po::value<bool>()->default_value(true),
      "Record the material interaction and - in extrapolation mode only.")(
      "prop-material-collection",
      po::value<std::string>()->default_value("propagation-material"),
      "Propagation material collection.")(
      "prop-ntests", po::value<std::size_t>()->default_value(1000),
      "Number of tests performed.")("prop-resolve-material",
                                    po::value<bool>()->default_value(true),
                                    "Resolve all smaterial surfaces.")(
      "prop-resolve-passive", po::value<bool>()->default_value(false),
      "Resolve all passive surfaces.")("prop-resolve-sensitive",
                                       po::value<bool>()->default_value(true),
                                       "Resolve all sensitive surfaces.")(
      "prop-d0-sigma", po::value<double>()->default_value(15_um),
      "Sigma of the transverse impact parameter [in mm].")(
      "prop-z0-sigma", po::value<double>()->default_value(55_mm),
      "Sigma of the longitudinal impact parameter [in mm].")(
      "prop-phi-sigma", po::value<double>()->default_value(0.001),
      "Sigma of the azimuthal angle [in rad].")(
      "prop-theta-sigma", po::value<double>()->default_value(0.001),
      "Sigma of the polar angle [in rad].")(
      "prop-qp-sigma", po::value<double>()->default_value(0.0001 / 1_GeV),
      "Sigma of the signed inverse momentum [in GeV^{-1}].")(
      "prop-t-sigma", po::value<double>()->default_value(1_ns),
      "Sigma of the time parameter [in ns].")(
      "prop-corr-offd", po::value<Reals<15>>(),
      "The 15 off-diagonal correlation rho(d0,z0), rho(d0,phi), [...], "
      "rho(z0,phi), rho(z0, theta), [...], rho(qop,t). Row-wise.")(
      "prop-phi-range", po::value<Reals<2>>()->default_value({{-M_PI, M_PI}}),
      "Azimutal angle phi range for proprapolated tracks.")(
      "prop-eta-range", po::value<Reals<2>>()->default_value({{-4., 4.}}),
      "Pseudorapidity range for proprapolated tracks.")(
      "prop-pt-range",
      po::value<Reals<2>>()->default_value({{100_MeV, 100_GeV}}),
      "Transverse momentum range for proprapolated tracks [in GeV].")(
      "prop-max-stepsize", po::value<double>()->default_value(3_m),
      "Maximum step size for the propagation [in mm].")(
      "prop-pt-loopers", po::value<double>()->default_value(500_MeV),
      "Transverse momentum below which loops are being detected [in GeV].");
}

/// Read the pgropagator options and return a Config file
///
/// @tparam vmap_t is the Type of the Parameter map to be read out
/// @tparam propagator_t is the Type of the Propagator used
///
/// @returns a Config object for the PropagationAlgorithm
inline ActsExamples::PropagationAlgorithm::Config readPropagationConfig(
    const boost::program_options::variables_map& vm) {
  using namespace Acts::UnitLiterals;
  ActsExamples::PropagationAlgorithm::Config pAlgConfig;

  auto iphir = vm["prop-phi-range"].template as<Reals<2>>();
  auto ietar = vm["prop-eta-range"].template as<Reals<2>>();
  auto iptr = vm["prop-pt-range"].template as<Reals<2>>();

  /// Material interaction behavior
  pAlgConfig.energyLoss = vm["prop-energyloss"].template as<bool>();
  pAlgConfig.multipleScattering = vm["prop-scattering"].template as<bool>();
  pAlgConfig.recordMaterialInteractions =
      vm["prop-record-material"].template as<bool>();

  /// Create the config for the Extrapoaltion algorithm
  pAlgConfig.debugOutput = vm["prop-debug"].template as<bool>();
  pAlgConfig.ntests = vm["prop-ntests"].template as<std::size_t>();
  pAlgConfig.mode = vm["prop-mode"].template as<int>();
  pAlgConfig.d0Sigma = vm["prop-d0-sigma"].template as<double>() * 1_mm;
  pAlgConfig.z0Sigma = vm["prop-z0-sigma"].template as<double>() * 1_mm;
  pAlgConfig.phiSigma = vm["prop-phi-sigma"].template as<double>();
  pAlgConfig.thetaSigma = vm["prop-theta-sigma"].template as<double>();
  pAlgConfig.qpSigma = vm["prop-qp-sigma"].template as<double>() / 1_GeV;
  pAlgConfig.tSigma = vm["prop-t-sigma"].template as<double>() * 1_ns;

  pAlgConfig.phiRange = {iphir[0], iphir[1]};
  pAlgConfig.etaRange = {ietar[0], ietar[1]};
  pAlgConfig.ptRange = {iptr[0] * 1_GeV, iptr[1] * 1_GeV};
  pAlgConfig.ptLoopers = vm["prop-pt-loopers"].template as<double>() * 1_GeV;
  pAlgConfig.maxStepSize = vm["prop-max-stepsize"].template as<double>() * 1_mm;

  pAlgConfig.propagationStepCollection =
      vm["prop-step-collection"].template as<std::string>();
  pAlgConfig.propagationMaterialCollection =
      vm["prop-material-collection"].template as<std::string>();

  /// The covariance transport
  if (vm["prop-cov"].template as<bool>()) {
    /// Set the covariance transport to true
    pAlgConfig.covarianceTransport = true;
    /// Set the covariance matrix
    pAlgConfig.covariances(Acts::BoundIndices::eBoundLoc0) =
        pAlgConfig.d0Sigma * pAlgConfig.d0Sigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundLoc1) =
        pAlgConfig.z0Sigma * pAlgConfig.z0Sigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundPhi) =
        pAlgConfig.phiSigma * pAlgConfig.phiSigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundTheta) =
        pAlgConfig.thetaSigma * pAlgConfig.thetaSigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundQOverP) =
        pAlgConfig.qpSigma * pAlgConfig.qpSigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundTime) =
        pAlgConfig.tSigma * pAlgConfig.tSigma;

    // Only if they are properly defined, assign off-diagonals
    if (vm.count("prop-corr-offd") != 0u) {
      auto readOffd = vm["prop-corr-offd"].template as<Reals<15>>();
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0,
                              Acts::BoundIndices::eBoundLoc1) = readOffd[0];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0,
                              Acts::BoundIndices::eBoundPhi) = readOffd[1];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0,
                              Acts::BoundIndices::eBoundTheta) = readOffd[2];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0,
                              Acts::BoundIndices::eBoundQOverP) = readOffd[3];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0,
                              Acts::BoundIndices::eBoundTime) = readOffd[4];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1,
                              Acts::BoundIndices::eBoundPhi) = readOffd[5];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1,
                              Acts::BoundIndices::eBoundTheta) = readOffd[6];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1,
                              Acts::BoundIndices::eBoundQOverP) = readOffd[7];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1,
                              Acts::BoundIndices::eBoundTime) = readOffd[8];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundPhi,
                              Acts::BoundIndices::eBoundTheta) = readOffd[9];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundPhi,
                              Acts::BoundIndices::eBoundQOverP) = readOffd[10];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundPhi,
                              Acts::BoundIndices::eBoundTime) = readOffd[11];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundTheta,
                              Acts::BoundIndices::eBoundQOverP) = readOffd[12];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundTheta,
                              Acts::BoundIndices::eBoundTime) = readOffd[13];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundQOverP,
                              Acts::BoundIndices::eBoundTime) = readOffd[14];
    } else {
      /// Some pre-defined values (non-trivial helical correlations)
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0,
                              Acts::BoundIndices::eBoundPhi) = -0.8;
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0,
                              Acts::BoundIndices::eBoundQOverP) = -0.3;
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1,
                              Acts::BoundIndices::eBoundTheta) = -0.8;
      pAlgConfig.correlations(Acts::BoundIndices::eBoundPhi,
                              Acts::BoundIndices::eBoundQOverP) = 0.4;
    }
  }

  return pAlgConfig;
}

}  // namespace Options
}  // namespace ActsExamples
