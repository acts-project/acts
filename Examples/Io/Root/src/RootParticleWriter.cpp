// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootParticleWriter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <cstdint>
#include <ios>
#include <ostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

ActsExamples::RootParticleWriter::RootParticleWriter(
    const ActsExamples::RootParticleWriter::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputParticles, "RootParticleWriter", lvl), m_cfg(cfg) {
  // inputParticles is already checked by base constructor
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // open root file and create the tree
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // setup the branches
  m_outputTree->Branch("event_id", &m_eventId);
  m_outputTree->Branch("particle_hash", &m_particleHash);
  m_outputTree->Branch("particle_type", &m_particleType);
  m_outputTree->Branch("process", &m_process);
  m_outputTree->Branch("vx", &m_vx);
  m_outputTree->Branch("vy", &m_vy);
  m_outputTree->Branch("vz", &m_vz);
  m_outputTree->Branch("vt", &m_vt);
  m_outputTree->Branch("px", &m_px);
  m_outputTree->Branch("py", &m_py);
  m_outputTree->Branch("pz", &m_pz);
  m_outputTree->Branch("m", &m_m);
  m_outputTree->Branch("q", &m_q);
  m_outputTree->Branch("eta", &m_eta);
  m_outputTree->Branch("phi", &m_phi);
  m_outputTree->Branch("pt", &m_pt);
  m_outputTree->Branch("p", &m_p);
  m_outputTree->Branch("q_over_p", &m_qop);
  m_outputTree->Branch("theta", &m_theta);
  m_outputTree->Branch("vertex_primary", &m_vertexPrimary);
  m_outputTree->Branch("vertex_secondary", &m_vertexSecondary);
  m_outputTree->Branch("particle", &m_particle);
  m_outputTree->Branch("generation", &m_generation);
  m_outputTree->Branch("sub_particle", &m_subParticle);

  if (m_cfg.writeHelixParameters) {
    m_outputTree->Branch("perigee_d0", &m_perigeeD0);
    m_outputTree->Branch("perigee_z0", &m_perigeeZ0);
    m_outputTree->Branch("perigee_phi", &m_perigeePhi);
    m_outputTree->Branch("perigee_theta", &m_perigeeTheta);
    m_outputTree->Branch("perigee_q_over_p", &m_perigeeQop);
    m_outputTree->Branch("perigee_p", &m_perigeeP);
    m_outputTree->Branch("perigee_px", &m_perigeePx);
    m_outputTree->Branch("perigee_py", &m_perigeePy);
    m_outputTree->Branch("perigee_pz", &m_perigeePz);
    m_outputTree->Branch("perigee_eta", &m_perigeeEta);
    m_outputTree->Branch("perigee_pt", &m_perigeePt);
  }

  m_outputTree->Branch("e_loss", &m_eLoss);
  m_outputTree->Branch("total_x0", &m_pathInX0);
  m_outputTree->Branch("total_l0", &m_pathInL0);
  m_outputTree->Branch("number_of_hits", &m_numberOfHits);
  m_outputTree->Branch("outcome", &m_outcome);
}

ActsExamples::RootParticleWriter::~RootParticleWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootParticleWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_INFO("Wrote particles to tree '" << m_cfg.treeName << "' in '"
                                        << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootParticleWriter::writeT(
    const AlgorithmContext& ctx, const SimParticleContainer& particles) {
  // ensure exclusive access to tree/file while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_eventId = ctx.eventNumber;
  for (const auto& particle : particles) {
    m_particleHash.push_back(particle.particleId().hash());
    m_particleType.push_back(particle.pdg());
    m_process.push_back(static_cast<std::uint32_t>(particle.process()));
    // position
    m_vx.push_back(Acts::clampValue<float>(particle.fourPosition().x() /
                                           Acts::UnitConstants::mm));
    m_vy.push_back(Acts::clampValue<float>(particle.fourPosition().y() /
                                           Acts::UnitConstants::mm));
    m_vz.push_back(Acts::clampValue<float>(particle.fourPosition().z() /
                                           Acts::UnitConstants::mm));
    m_vt.push_back(Acts::clampValue<float>(particle.fourPosition().w() /
                                           Acts::UnitConstants::mm));

    // particle constants
    if (!std::isfinite(particle.mass()) || !std::isfinite(particle.charge())) {
      ACTS_WARNING("Particle mass or charge is not finite, can't write it");
    }

    m_m.push_back(
        Acts::clampValue<float>(particle.mass() / Acts::UnitConstants::GeV));
    m_q.push_back(
        Acts::clampValue<float>(particle.charge() / Acts::UnitConstants::e));
    // decoded barcode components
    m_vertexPrimary.push_back(particle.particleId().vertexPrimary());
    m_vertexSecondary.push_back(particle.particleId().vertexSecondary());
    m_particle.push_back(particle.particleId().particle());
    m_generation.push_back(particle.particleId().generation());
    m_subParticle.push_back(particle.particleId().subParticle());

    m_eLoss.push_back(Acts::clampValue<float>(particle.energyLoss() /
                                              Acts::UnitConstants::GeV));
    m_pathInX0.push_back(
        Acts::clampValue<float>(particle.pathInX0() / Acts::UnitConstants::mm));
    m_pathInL0.push_back(
        Acts::clampValue<float>(particle.pathInL0() / Acts::UnitConstants::mm));
    m_numberOfHits.push_back(particle.numberOfHits());
    m_outcome.push_back(static_cast<std::uint32_t>(particle.outcome()));

    // momentum
    const auto p = particle.absoluteMomentum() / Acts::UnitConstants::GeV;
    m_p.push_back(Acts::clampValue<float>(p));
    m_px.push_back(Acts::clampValue<float>(p * particle.direction().x()));
    m_py.push_back(Acts::clampValue<float>(p * particle.direction().y()));
    m_pz.push_back(Acts::clampValue<float>(p * particle.direction().z()));
    // derived kinematic quantities
    m_eta.push_back(Acts::clampValue<float>(
        Acts::VectorHelpers::eta(particle.direction())));
    m_pt.push_back(Acts::clampValue<float>(
        p * Acts::VectorHelpers::perp(particle.direction())));
    m_phi.push_back(Acts::clampValue<float>(
        Acts::VectorHelpers::phi(particle.direction())));
    m_theta.push_back(Acts::clampValue<float>(
        Acts::VectorHelpers::theta(particle.direction())));
    m_qop.push_back(Acts::clampValue<float>(
        particle.qOverP() * Acts::UnitConstants::GeV / Acts::UnitConstants::e));

    if (!m_cfg.writeHelixParameters) {
      // done with this particle
      continue;
    }

    // Perigee surface at configured reference point
    auto pSurface =
        Acts::Surface::makeShared<Acts::PerigeeSurface>(m_cfg.referencePoint);

    // Start from truth curvilinear parameters (direction, q/p)
    const Acts::Vector3 startDir = particle.direction();  // unit vector
    const auto qOverP = particle.qOverP();                // ACTS units
    auto intersection =
        pSurface
            ->intersect(ctx.geoContext, particle.position(), startDir,
                        Acts::BoundaryTolerance::Infinite())
            .closest();

    // Neutral particles have no helix -> linearly extrapolate to perigee
    if (particle.charge() == 0) {
      ACTS_WARNING(
          "Particle has zero charge, linearly extrapolating to perigee");
      // Initialize the truth particle info
      auto perigeeD0 = NaNfloat;
      auto perigeeZ0 = NaNfloat;

      const auto position = intersection.position();

      // get the truth perigee parameter
      auto lpResult =
          pSurface->globalToLocal(ctx.geoContext, position, startDir);
      if (lpResult.ok()) {
        perigeeD0 = lpResult.value()[Acts::BoundIndices::eBoundLoc0];
        perigeeZ0 = lpResult.value()[Acts::BoundIndices::eBoundLoc1];
      } else {
        ACTS_ERROR("Global to local transformation did not succeed.");
      }
      // truth parameters at perigee are the same as at production vertex
      m_perigeePhi.push_back(Acts::clampValue<float>(particle.phi()));
      m_perigeeTheta.push_back(Acts::clampValue<float>(particle.theta()));
      m_perigeeQop.push_back(Acts::clampValue<float>(
          qOverP * Acts::UnitConstants::GeV / Acts::UnitConstants::e));
      m_perigeeP.push_back(Acts::clampValue<float>(particle.absoluteMomentum() /
                                                   Acts::UnitConstants::GeV));
      m_perigeePx.push_back(Acts::clampValue<float>(m_p.back() * startDir.x()));
      m_perigeePy.push_back(Acts::clampValue<float>(m_p.back() * startDir.y()));
      m_perigeePz.push_back(Acts::clampValue<float>(m_p.back() * startDir.z()));
      m_perigeeEta.push_back(Acts::clampValue<float>(
          Acts::VectorHelpers::eta(particle.direction())));
      m_perigeePt.push_back(Acts::clampValue<float>(
          m_p.back() * Acts::VectorHelpers::perp(particle.direction())));

      // Push the extrapolated parameters
      m_perigeeD0.push_back(
          Acts::clampValue<float>(perigeeD0 / Acts::UnitConstants::mm));
      m_perigeeZ0.push_back(
          Acts::clampValue<float>(perigeeZ0 / Acts::UnitConstants::mm));
      continue;
    }

    // Charged particles: propagate helix to perigee
    // Build a propagator and propagate the *truth parameters* to the
    // perigee Stepper + propagator
    using Stepper = Acts::SympyStepper;
    Stepper stepper(m_cfg.bField);
    using PropagatorT = Acts::Propagator<Stepper>;
    auto propagator = std::make_shared<PropagatorT>(stepper);

    Acts::BoundTrackParameters startParams =
        Acts::BoundTrackParameters::createCurvilinear(
            particle.fourPosition(), startDir, qOverP, std::nullopt,
            Acts::ParticleHypothesis::pion());

    // Propagation options (need event contexts)
    using PropOptions = PropagatorT::Options<>;
    PropOptions pOptions(ctx.geoContext, ctx.magFieldContext);

    // Choose propagation direction based on the closest intersection
    pOptions.direction =
        Acts::Direction::fromScalarZeroAsPositive(intersection.pathLength());

    // Do the propagation to the perigee surface
    auto propRes = propagator->propagate(startParams, *pSurface, pOptions);
    if (!propRes.ok() || !propRes->endParameters.has_value()) {
      ACTS_ERROR("Propagation to perigee surface failed.");
      m_perigeePhi.push_back(NaNfloat);
      m_perigeeTheta.push_back(NaNfloat);
      m_perigeeQop.push_back(NaNfloat);
      m_perigeeD0.push_back(NaNfloat);
      m_perigeeZ0.push_back(NaNfloat);
      m_perigeeP.push_back(NaNfloat);
      m_perigeePx.push_back(NaNfloat);
      m_perigeePy.push_back(NaNfloat);
      m_perigeePz.push_back(NaNfloat);
      m_perigeeEta.push_back(NaNfloat);
      m_perigeePt.push_back(NaNfloat);
      continue;
    }
    const Acts::BoundTrackParameters& atPerigee = *propRes->endParameters;

    // By construction, atPerigee is *bound on the perigee surface*.
    // Its parameter vector is [loc0, loc1, phi, theta, q/p, t]
    const auto& perigee_pars = atPerigee.parameters();

    const auto perigeeD0 = perigee_pars[Acts::BoundIndices::eBoundLoc0];
    const auto perigeeZ0 = perigee_pars[Acts::BoundIndices::eBoundLoc1];

    // truth phi;theta;q/p *at the perigee*,
    const auto perigeePhi = perigee_pars[Acts::BoundIndices::eBoundPhi];
    const auto perigeeTheta = perigee_pars[Acts::BoundIndices::eBoundTheta];
    const auto perigeeQop = perigee_pars[Acts::BoundIndices::eBoundQOverP];

    m_perigeePhi.push_back(Acts::clampValue<float>(perigeePhi));
    m_perigeeTheta.push_back(Acts::clampValue<float>(perigeeTheta));
    m_perigeeQop.push_back(Acts::clampValue<float>(
        perigeeQop * Acts::UnitConstants::GeV / Acts::UnitConstants::e));
    // update p, px, py, pz, eta, pt
    const auto perigeeP =
        atPerigee.absoluteMomentum() / Acts::UnitConstants::GeV;
    m_perigeeP.push_back(Acts::clampValue<float>(perigeeP));
    const auto dir = atPerigee.direction();
    m_perigeePx.push_back(Acts::clampValue<float>(perigeeP * dir.x()));
    m_perigeePy.push_back(Acts::clampValue<float>(perigeeP * dir.y()));
    m_perigeePz.push_back(Acts::clampValue<float>(perigeeP * dir.z()));
    m_perigeeEta.push_back(Acts::clampValue<float>(
        Acts::VectorHelpers::eta(atPerigee.direction())));
    m_perigeePt.push_back(Acts::clampValue<float>(
        perigeeP * Acts::VectorHelpers::perp(atPerigee.direction())));

    // Push the perigee parameters
    m_perigeeD0.push_back(
        Acts::clampValue<float>(perigeeD0 / Acts::UnitConstants::mm));
    m_perigeeZ0.push_back(
        Acts::clampValue<float>(perigeeZ0 / Acts::UnitConstants::mm));
  }

  m_outputTree->Fill();

  m_particleHash.clear();
  m_particleType.clear();
  m_process.clear();
  m_vx.clear();
  m_vy.clear();
  m_vz.clear();
  m_vt.clear();
  m_p.clear();
  m_px.clear();
  m_py.clear();
  m_pz.clear();
  m_m.clear();
  m_q.clear();
  m_eta.clear();
  m_phi.clear();
  m_pt.clear();
  m_theta.clear();
  m_qop.clear();
  m_vertexPrimary.clear();
  m_vertexSecondary.clear();
  m_particle.clear();
  m_generation.clear();
  m_subParticle.clear();
  m_eLoss.clear();
  m_numberOfHits.clear();
  m_outcome.clear();
  m_pathInX0.clear();
  m_pathInL0.clear();

  if (m_cfg.writeHelixParameters) {
    m_perigeeD0.clear();
    m_perigeeZ0.clear();
    m_perigeePhi.clear();
    m_perigeeTheta.clear();
    m_perigeeQop.clear();
    m_perigeeP.clear();
    m_perigeePx.clear();
    m_perigeePy.clear();
    m_perigeePz.clear();
    m_perigeeEta.clear();
    m_perigeePt.clear();
  }

  return ProcessCode::SUCCESS;
}
