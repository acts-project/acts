// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"

#include <stdexcept>

namespace ActsExamples {

ProcessCode PropagationAlgorithm::execute(
    const AlgorithmContext& context) const {
  // Create a random number generator
  ActsExamples::RandomEngine rng =
      m_cfg.randomNumberSvc->spawnGenerator(context);

  // Standard gaussian distribution for covarianmces
  std::normal_distribution<double> gauss(0., 1.);

  // Setup random number distributions for some quantities
  std::uniform_real_distribution<double> phiDist(m_cfg.phiRange.first,
                                                 m_cfg.phiRange.second);
  std::uniform_real_distribution<double> etaDist(m_cfg.etaRange.first,
                                                 m_cfg.etaRange.second);
  std::uniform_real_distribution<double> ptDist(m_cfg.ptRange.first,
                                                m_cfg.ptRange.second);
  std::uniform_real_distribution<double> qDist(0., 1.);

  std::shared_ptr<const Acts::PerigeeSurface> surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(
          Acts::Vector3(0., 0., 0.));

  // Output : the propagation steps
  std::vector<std::vector<Acts::detail::Step>> propagationSteps;
  propagationSteps.reserve(m_cfg.ntests);

  // Output (optional): the recorded material
  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack> recordedMaterial;

  // loop over number of particles
  for (std::size_t it = 0; it < m_cfg.ntests; ++it) {
    /// get the d0 and z0
    double d0 = m_cfg.d0Sigma * gauss(rng);
    double z0 = m_cfg.z0Sigma * gauss(rng);
    double phi = phiDist(rng);
    double eta = etaDist(rng);
    double theta = 2 * std::atan(std::exp(-eta));
    double pt = ptDist(rng);
    double p = pt / std::sin(theta);
    double charge = qDist(rng) > 0.5 ? 1. : -1.;
    double qop = charge / p;
    double t = m_cfg.tSigma * gauss(rng);
    // parameters
    Acts::BoundVector pars;
    pars << d0, z0, phi, theta, qop, t;
    // some screen output

    // The covariance generation
    auto cov = generateCovariance(rng, gauss);

    // execute the test for charged particles
    Acts::BoundTrackParameters startParameters(surface, pars, std::move(cov),
                                               m_cfg.particleHypothesis);
    Acts::Vector3 sPosition = startParameters.position(context.geoContext);
    Acts::Vector3 sMomentum = startParameters.momentum();
    PropagationOutput pOutput = m_cfg.propagatorImpl->execute(
        context, m_cfg, logger(), startParameters);
    // Record the propagator steps
    propagationSteps.push_back(std::move(pOutput.first));
    if (m_cfg.recordMaterialInteractions &&
        !pOutput.second.materialInteractions.empty()) {
      // Create a recorded material track with start position, momentum and the
      // material
      recordedMaterial.emplace(
          it, std::make_pair(std::make_pair(sPosition, sMomentum),
                             std::move(pOutput.second)));
    }
  }

  // Write the propagation step data to the event store
  m_outpoutPropagationSteps(context, std::move(propagationSteps));

  // Write the recorded material to the event store
  if (m_cfg.recordMaterialInteractions) {
    m_recordedMaterial(context, std::move(recordedMaterial));
  }

  return ProcessCode::SUCCESS;
}

std::optional<Acts::BoundSquareMatrix> PropagationAlgorithm::generateCovariance(
    ActsExamples::RandomEngine& rnd,
    std::normal_distribution<double>& gauss) const {
  if (m_cfg.covarianceTransport) {
    // We start from the correlation matrix
    Acts::BoundSquareMatrix newCov(m_cfg.correlations);
    // Then we draw errors according to the error values
    Acts::BoundVector covs_smeared = m_cfg.covariances;
    for (std::size_t k = 0; k < static_cast<std::size_t>(covs_smeared.size());
         ++k) {
      covs_smeared[k] *= gauss(rnd);
    }
    // and apply a double loop
    for (std::size_t i = 0; i < static_cast<std::size_t>(newCov.rows()); ++i) {
      for (std::size_t j = 0; j < static_cast<std::size_t>(newCov.cols());
           ++j) {
        (newCov)(i, j) *= covs_smeared[i];
        (newCov)(i, j) *= covs_smeared[j];
      }
    }
    return newCov;
  }
  return std::nullopt;
}

PropagationAlgorithm::PropagationAlgorithm(
    const PropagationAlgorithm::Config& config, Acts::Logging::Level level)
    : IAlgorithm("PropagationAlgorithm", level), m_cfg(config) {
  if (!m_cfg.propagatorImpl) {
    throw std::invalid_argument("Config needs to contain a propagator");
  }
  if (!m_cfg.randomNumberSvc) {
    throw std::invalid_argument("No random number generator given");
  }

  m_outpoutPropagationSteps.initialize(m_cfg.propagationStepCollection);
  m_recordedMaterial.initialize(m_cfg.propagationMaterialCollection);
}

}  // namespace ActsExamples
