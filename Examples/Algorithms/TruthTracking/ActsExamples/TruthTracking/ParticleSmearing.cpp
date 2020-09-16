// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <cmath>
#include <vector>

ActsExamples::ParticleSmearing::ParticleSmearing(const Config& cfg,
                                                 Acts::Logging::Level lvl)
    : BareAlgorithm("ParticleSmearing", lvl), m_cfg(cfg) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing output tracks parameters collection");
  }
}

ActsExamples::ProcessCode ActsExamples::ParticleSmearing::execute(
    const AlgorithmContext& ctx) const {
  // setup input and output containers
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  TrackParametersContainer parameters;
  parameters.reserve(particles.size());

  // setup random number generator and standard gaussian
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  std::normal_distribution<double> stdNormal(0.0, 1.0);

  for (const auto& particle : particles) {
    const auto phi = Acts::VectorHelpers::phi(particle.unitDirection());
    const auto theta = Acts::VectorHelpers::theta(particle.unitDirection());
    const auto pt = particle.transverseMomentum();
    const auto p = particle.absMomentum();

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

    // smear the position/time
    Acts::Vector4D pos4 = particle.position4();
    pos4[Acts::ePos0] += deltaD0 * std::sin(phi);
    pos4[Acts::ePos1] += deltaD0 * -std::cos(phi);
    pos4[Acts::ePos2] += deltaZ0;
    pos4[Acts::eTime] += deltaT0;
    // smear direction angles phi,theta ensuring correct bounds
    const auto [newPhi, newTheta] =
        Acts::detail::ensureThetaBounds(phi + deltaPhi, theta + deltaTheta);
    // compute smeared absolute momentum vector
    const auto newP = p + deltaP;

    // build the track covariance matrix using the smearing sigmas
    Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = sigmaU * sigmaU;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = sigmaV * sigmaV;
    cov(Acts::eBoundTime, Acts::eBoundTime) = sigmaT0 * sigmaT0;
    cov(Acts::eBoundPhi, Acts::eBoundPhi) = sigmaPhi * sigmaPhi;
    cov(Acts::eBoundTheta, Acts::eBoundTheta) = sigmaTheta * sigmaTheta;
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = sigmaQOverP * sigmaQOverP;

    parameters.emplace_back(pos4, newPhi, newTheta, newP, particle.charge(),
                            cov);
  };

  ctx.eventStore.add(m_cfg.outputTrackParameters, std::move(parameters));
  return ProcessCode::SUCCESS;
}
