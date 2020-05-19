// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/TruthTracking/TruthVerticesToTracks.hpp"

#include <iostream>
#include <optional>
#include <stdexcept>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"

FW::TruthVerticesToTracksAlgorithm::TruthVerticesToTracksAlgorithm(
    const FW::TruthVerticesToTracksAlgorithm::Config& cfg,
    Acts::Logging::Level level)
    : FW::BareAlgorithm("TruthVerticesToTracksAlgorithm", level), m_cfg(cfg) {
  if (m_cfg.input.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.output.empty()) {
    throw std::invalid_argument("Missing output collection");
  } else if (m_cfg.randomNumberSvc == nullptr) {
    throw std::invalid_argument("Missing random number service");
  }
}

FW::ProcessCode FW::TruthVerticesToTracksAlgorithm::execute(
    const AlgorithmContext& context) const {
  const auto& vertexCollection =
      context.eventStore.get<std::vector<FW::SimVertex>>(m_cfg.input);

  std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(m_cfg.refPosition);

  // Set up constant B-Field
  Acts::ConstantBField bField(m_cfg.bField);

  // Set up stepper
  Acts::EigenStepper<Acts::ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  Acts::Propagator<Acts::EigenStepper<Acts::ConstantBField>> propagator(
      stepper);

  // Set up propagator options
  Acts::PropagatorOptions<> pOptions(context.geoContext,
                                     context.magFieldContext);
  pOptions.direction = Acts::backward;

  // Create random number generator and spawn gaussian distribution
  FW::RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);

  // Vector to store VertexAndTracks extracted from event
  std::vector<VertexAndTracks> vertexAndTracksCollection;

  // Start looping over all vertices in current event
  for (auto& vtx : vertexCollection) {
    // Create VertexAndTracks object
    VertexAndTracks vertexAndTracks;
    // Store current vertex
    vertexAndTracks.vertex = vtx;

    // Track objects at current vertex
    std::vector<Acts::BoundParameters> trackCollection;

    // Iterate over all particle emerging from current vertex
    for (auto const& particle : vtx.outgoing) {
      const Acts::Vector3D& ptclMom =
          particle.absMomentum() * particle.unitDirection();

      // Define start track params
      Acts::CurvilinearParameters start(std::nullopt, particle.position(),
                                        ptclMom, particle.charge(),
                                        particle.time());
      // Run propagator
      auto result = propagator.propagate(start, *perigeeSurface, pOptions);
      if (!result.ok()) {
        continue;
      }

      // get perigee parameters
      const auto& perigeeParameters = (*result).endParameters->parameters();

      auto newTrackParams = perigeeParameters;

      if (m_cfg.doSmearing) {
        // Calculate pt-dependent IP resolution
        const double particlePt = Acts::VectorHelpers::perp(ptclMom);
        const double ipRes =
            m_cfg.ipResA * std::exp(-m_cfg.ipResB * particlePt) + m_cfg.ipResC;

        // except for IP resolution, following variances are rough guesses
        // Gaussian distribution for IP resolution
        std::normal_distribution<double> gaussDist_IP(0., ipRes);
        // Gaussian distribution for angular resolution
        std::normal_distribution<double> gaussDist_angular(0., m_cfg.angRes);
        // Gaussian distribution for q/p (momentum) resolution
        std::normal_distribution<double> gaussDist_qp(
            0., m_cfg.qpRelRes * perigeeParameters[4]);

        double rn_d0 = gaussDist_IP(rng);
        double rn_z0 = gaussDist_IP(rng);
        double rn_ph = gaussDist_angular(rng);
        double rn_th = gaussDist_angular(rng);
        double rn_qp = gaussDist_qp(rng);

        Acts::BoundVector smrdParamVec;
        smrdParamVec << rn_d0, rn_z0, rn_ph, rn_th, rn_qp, 0.;

        // Update track parameters
        newTrackParams += smrdParamVec;

        // Correct for phi and theta wrap
        correctPhiThetaPeriodicity(newTrackParams[2], newTrackParams[3]);

        // Update track covariance
        Acts::BoundSymMatrix covMat;
        covMat.setZero();
        covMat.diagonal() << rn_d0 * rn_d0, rn_z0 * rn_z0, rn_ph * rn_ph,
            rn_th * rn_th, rn_qp * rn_qp, 1.;

        trackCollection.push_back(Acts::BoundParameters(
            context.geoContext, covMat, newTrackParams, perigeeSurface));
      } else {
        trackCollection.push_back(Acts::BoundParameters(
            context.geoContext, std::nullopt, newTrackParams, perigeeSurface));
      }
    }  // end iteration over all particle at vertex

    // Store track objects in VertexAndTracks
    vertexAndTracks.tracks = trackCollection;
    // Add to collection
    vertexAndTracksCollection.push_back(vertexAndTracks);

  }  // end iteration over all vertices

  // VertexAndTracks objects to the EventStore
  context.eventStore.add(m_cfg.output, std::move(vertexAndTracksCollection));

  return FW::ProcessCode::SUCCESS;
}

void FW::TruthVerticesToTracksAlgorithm::correctPhiThetaPeriodicity(
    double& phiIn, double& thetaIn) const {
  double tmpPhi = std::fmod(phiIn, 2 * M_PI);  // temp phi
  if (tmpPhi > M_PI) {
    tmpPhi -= 2 * M_PI;
  }
  if (tmpPhi < -M_PI && tmpPhi > -2 * M_PI) {
    tmpPhi += 2 * M_PI;
  }

  double tmpTht = std::fmod(thetaIn, 2 * M_PI);  // temp theta
  if (tmpTht < -M_PI) {
    tmpTht = std::abs(tmpTht + 2 * M_PI);
  } else if (tmpTht < 0) {
    tmpTht *= -1;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? tmpPhi - 2 * M_PI : tmpPhi;
  }
  if (tmpTht > M_PI) {
    tmpTht = 2 * M_PI - tmpTht;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? (tmpPhi - 2 * M_PI) : tmpPhi;
  }

  phiIn = tmpPhi;
  thetaIn = tmpTht;
}
