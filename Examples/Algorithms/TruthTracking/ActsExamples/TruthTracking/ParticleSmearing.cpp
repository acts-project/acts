// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <cmath>
#include <vector>

#include <boost/program_options.hpp>

void ActsExamples::ParticleSmearing::addOptions(Options::Description& desc) {
  using boost::program_options::value;
  using Options::Reals;

  auto opt = desc.add_options();
  opt("smear-sigma-D0", value<Reals<3>>()->default_value({{20, 30, 0.3}}),
      "Smear the initial Pt-dependent d0 in perigee frame with a_0[um] + "
      "a_1[um]*exp(-1.*abs(a_2[1/GeV])*pt)");
  opt("smear-sigma-Z0", value<Reals<3>>()->default_value({{20, 30, 0.3}}),
      "Smear the initial Pt-dependent z0 in perigee frame with a_0[um] + "
      "a_1[um]*exp(-1.*abs(a_2[1/GeV])*pt)");
  opt("smear-sigma-T0", value<double>()->default_value(1),
      "Smear the initial time in ns");
  opt("smear-sigma-momentum", value<Reals<3>>()->default_value({{1, 1, 0.1}}),
      "Smear the initial phi (degree), theta (degree) and momentum (relative)");
}

ActsExamples::ParticleSmearing::Config
ActsExamples::ParticleSmearing::readConfig(const Options::Variables& vars) {
  using namespace Acts::UnitConstants;
  using Options::Reals;

  Config cfg;
  auto sigmaD0 = vars["smear-sigma-D0"].template as<Reals<3>>();
  auto sigmaZ0 = vars["smear-sigma-Z0"].template as<Reals<3>>();
  auto sigmaMom = vars["smear-sigma-momentum"].template as<Reals<3>>();

  cfg.sigmaD0 = sigmaD0[0] * Acts::UnitConstants::um;
  cfg.sigmaD0PtA = sigmaD0[1] * Acts::UnitConstants::um;
  cfg.sigmaD0PtB = sigmaD0[2] / Acts::UnitConstants::GeV;
  cfg.sigmaZ0 = sigmaZ0[0] * Acts::UnitConstants::um;
  cfg.sigmaZ0PtA = sigmaZ0[1] * Acts::UnitConstants::um;
  cfg.sigmaZ0PtB = sigmaZ0[2] / Acts::UnitConstants::GeV;
  cfg.sigmaT0 = vars["smear-sigma-T0"].as<double>() * Acts::UnitConstants::ns;
  cfg.sigmaPhi = sigmaMom[0] * Acts::UnitConstants::degree;
  cfg.sigmaTheta = sigmaMom[1] * Acts::UnitConstants::degree;
  cfg.sigmaPRel = sigmaMom[2];

  return cfg;
}

ActsExamples::ParticleSmearing::ParticleSmearing(const Config& config,
                                                 Acts::Logging::Level level)
    : BareAlgorithm("ParticleSmearing", level), m_cfg(config) {
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

  for (auto&& [vtxId, vtxParticles] : groupBySecondaryVertex(particles)) {
    // a group contains at least one particle by construction. assume that all
    // particles within the group originate from the same position and use it to
    // as the refernce position for the perigee frame.
    auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        vtxParticles.begin()->position());

    for (const auto& particle : vtxParticles) {
      const auto time = particle.time();
      const auto phi = Acts::VectorHelpers::phi(particle.unitDirection());
      const auto theta = Acts::VectorHelpers::theta(particle.unitDirection());
      const auto pt = particle.transverseMomentum();
      const auto p = particle.absoluteMomentum();
      const auto q = particle.charge();

      // compute momentum-dependent resolutions
      const double sigmaD0 =
          m_cfg.sigmaD0 +
          m_cfg.sigmaD0PtA * std::exp(-1.0 * std::abs(m_cfg.sigmaD0PtB) * pt);
      const double sigmaZ0 =
          m_cfg.sigmaZ0 +
          m_cfg.sigmaZ0PtA * std::exp(-1.0 * std::abs(m_cfg.sigmaZ0PtB) * pt);
      const double sigmaP = m_cfg.sigmaPRel * p;
      // var(q/p) = (d(1/p)/dp)² * var(p) = (-1/p²)² * var(p)
      const double sigmaQOverP = sigmaP / (p * p);
      // shortcuts for other resolutions
      const double sigmaT0 = m_cfg.sigmaT0;
      const double sigmaPhi = m_cfg.sigmaPhi;
      const double sigmaTheta = m_cfg.sigmaTheta;

      Acts::BoundVector params = Acts::BoundVector::Zero();
      // smear the position/time
      params[Acts::eBoundLoc0] = sigmaD0 * stdNormal(rng);
      params[Acts::eBoundLoc1] = sigmaZ0 * stdNormal(rng);
      params[Acts::eBoundTime] = time + sigmaT0 * stdNormal(rng);
      // smear direction angles phi,theta ensuring correct bounds
      const auto [newPhi, newTheta] = Acts::detail::normalizePhiTheta(
          phi + sigmaPhi * stdNormal(rng), theta + sigmaTheta * stdNormal(rng));
      params[Acts::eBoundPhi] = newPhi;
      params[Acts::eBoundTheta] = newTheta;
      // compute smeared absolute momentum vector
      const double newP = std::max(0.0, p + sigmaP * stdNormal(rng));
      params[Acts::eBoundQOverP] = (q != 0) ? (q / newP) : (1 / newP);

      ACTS_VERBOSE("Smearing particle (pos, time, phi, theta, q/p):");
      ACTS_VERBOSE(" from: " << particle.position().transpose() << ", " << time
                             << "," << phi << "," << theta << ","
                             << (q != 0 ? q / p : 1 / p));
      ACTS_VERBOSE("   to: " << perigee
                                    ->localToGlobal(
                                        ctx.geoContext,
                                        Acts::Vector2{params[Acts::eBoundLoc0],
                                                      params[Acts::eBoundLoc1]},
                                        particle.unitDirection() * p)
                                    .transpose()
                             << ", " << params[Acts::eBoundTime] << ","
                             << params[Acts::eBoundPhi] << ","
                             << params[Acts::eBoundTheta] << ","
                             << params[Acts::eBoundQOverP]);

      // build the track covariance matrix using the smearing sigmas
      Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) =
          m_cfg.initialVarInflation[Acts::eBoundLoc0] * sigmaD0 * sigmaD0;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) =
          m_cfg.initialVarInflation[Acts::eBoundLoc1] * sigmaZ0 * sigmaZ0;
      cov(Acts::eBoundTime, Acts::eBoundTime) =
          m_cfg.initialVarInflation[Acts::eBoundTime] * sigmaT0 * sigmaT0;
      cov(Acts::eBoundPhi, Acts::eBoundPhi) =
          m_cfg.initialVarInflation[Acts::eBoundPhi] * sigmaPhi * sigmaPhi;
      cov(Acts::eBoundTheta, Acts::eBoundTheta) =
          m_cfg.initialVarInflation[Acts::eBoundTheta] * sigmaTheta *
          sigmaTheta;
      cov(Acts::eBoundQOverP, Acts::eBoundQOverP) =
          m_cfg.initialVarInflation[Acts::eBoundQOverP] * sigmaQOverP *
          sigmaQOverP;

      parameters.emplace_back(perigee, params, q, cov);
    }
  }

  ctx.eventStore.add(m_cfg.outputTrackParameters, std::move(parameters));
  return ProcessCode::SUCCESS;
}
