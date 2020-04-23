// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/TruthTracking/ParticleSmearing.hpp"

#include <cmath>
#include <vector>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

FW::ParticleSmearing::ParticleSmearing(const Config& cfg,
                                       Acts::Logging::Level lvl)
    : BareAlgorithm("ParticleSmearing", lvl), m_cfg(cfg) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing output tracks parameters collection");
  }
}

FW::ProcessCode FW::ParticleSmearing::execute(
    const AlgorithmContext& ctx) const {
  namespace vh = Acts::VectorHelpers;

  // setup input and output containers
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  TrackParametersContainer parameters;
  parameters.reserve(particles.size());

  // setup random number generator and standard gaussian
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  std::normal_distribution<double> stdNormal(0.0, 1.0);

  for (const auto& particle : particles) {
    const auto pt = particle.transverseMomentum();
    const auto p = particle.absMomentum();
    const auto theta = vh::theta(particle.unitDirection());
    const auto phi = vh::phi(particle.unitDirection());

    // compute momentum-dependent resolutions
    const auto sigmaD0 =
        m_cfg.sigmaD0 +
        m_cfg.sigmaD0PtA * std::exp(-1.0 * std::abs(m_cfg.sigmaD0PtB) * pt);
    const auto sigmaZ0 =
        m_cfg.sigmaZ0 +
        m_cfg.sigmaZ0PtA * std::exp(-1.0 * std::abs(m_cfg.sigmaZ0PtB) * pt);
    const auto sigmaP = m_cfg.sigmaPRel * p;
    // var(q/p) = (d(1/p)/dp)² * var(p) = (-1/p²)² * var(p)
    const auto sigmaQOverP = sigmaP / (p * p);
    // shortcuts for other resolutions
    const auto sigmaT0 = m_cfg.sigmaT0;
    const auto sigmaPhi = m_cfg.sigmaPhi;
    const auto sigmaTheta = m_cfg.sigmaTheta;
    // converstion from perigee d0,z0 to curvilinear u,v
    // d0 and u differ only by a sign
    const auto sigmaU = sigmaD0;
    // project from z0 to the second axes orthogonal to the track direction
    const auto sigmaV = sigmaZ0 * std::sin(theta);

    // draw random noise
    const auto deltaD0 = sigmaD0 * stdNormal(rng);
    const auto deltaZ0 = sigmaZ0 * stdNormal(rng);
    const auto deltaT0 = sigmaT0 * stdNormal(rng);
    const auto deltaPhi = sigmaPhi * stdNormal(rng);
    const auto deltaTheta = sigmaTheta * stdNormal(rng);
    const auto deltaP = sigmaP * stdNormal(rng);

    // smear the position
    const Acts::Vector3D pos =
        particle.position() + Acts::Vector3D(deltaD0 * std::sin(phi),
                                             deltaD0 * -std::cos(phi), deltaZ0);
    // smear the time
    const auto time = particle.time() + deltaT0;
    // smear direction angles phi,theta ensuring correct bounds
    const auto angles =
        Acts::detail::ensureThetaBounds(phi + deltaPhi, theta + deltaTheta);
    // compute smeared direction vector
    const Acts::Vector3D dir(std::sin(angles.second) * std::cos(angles.first),
                             std::sin(angles.second) * std::sin(angles.first),
                             std::cos(angles.second));
    // compute smeared momentum vector
    const Acts::Vector3D mom = (p + deltaP) * dir;

    // build the track covariance matrix using the smearing sigmas
    Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
    cov(Acts::eLOC_0, Acts::eLOC_0) = sigmaU * sigmaU;
    cov(Acts::eLOC_1, Acts::eLOC_1) = sigmaV * sigmaV;
    cov(Acts::ePHI, Acts::ePHI) = sigmaPhi * sigmaPhi;
    cov(Acts::eTHETA, Acts::eTHETA) = sigmaTheta * sigmaTheta;
    cov(Acts::eQOP, Acts::eQOP) = sigmaQOverP * sigmaQOverP;
    cov(Acts::eT, Acts::eT) = sigmaT0 * sigmaT0;

    parameters.emplace_back(std::make_optional(std::move(cov)), pos, mom,
                            particle.charge(), time);
  };

  ctx.eventStore.add(m_cfg.outputTrackParameters, std::move(parameters));
  return ProcessCode::SUCCESS;
}
