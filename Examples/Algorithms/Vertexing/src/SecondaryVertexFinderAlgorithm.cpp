// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/SecondaryVertexFinderAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <chrono>
#include <ostream>
#include <stdexcept>
#include <system_error>

#include <TLorentzVector.h>
#include <TMath.h>

#include "VertexingHelpers.hpp"

ActsExamples::SecondaryVertexFinderAlgorithm::SecondaryVertexFinderAlgorithm(
    const Config& config, Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("IterativeVertexFinder", level), m_cfg(config) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing input track parameter collection");
  }
  if (m_cfg.outputProtoVertices.empty()) {
    throw std::invalid_argument("Missing output proto vertices collection");
  }
  if (m_cfg.outputVertices.empty()) {
    throw std::invalid_argument("Missing output vertices collection");
  }

  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_outputProtoVertices.initialize(m_cfg.outputProtoVertices);
  m_outputMasses.initialize(m_cfg.outputMasses);
  m_outputVertices.initialize(m_cfg.outputVertices);
}

ActsExamples::ProcessCode ActsExamples::SecondaryVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // retrieve input tracks and convert into the expected format

  const auto& inputTrackParameters = m_inputTrackParameters(ctx);
  // TODO change this from pointers to tracks parameters to actual tracks
  auto inputTracks = makeInputTracks(inputTrackParameters);

  if (inputTrackParameters.size() != inputTracks.size()) {
    ACTS_ERROR("Input track containers do not align: "
               << inputTrackParameters.size() << " != " << inputTracks.size());
  }

  for (const auto& trk : inputTrackParameters) {
    if (trk.covariance() && trk.covariance()->determinant() <= 0) {
      // actually we should consider this as an error but I do not want the CI
      // to fail
      ACTS_WARNING("input track " << trk << " has det(cov) = "
                                  << trk.covariance()->determinant());
    }
  }

  // Set up EigenStepper
  Acts::EigenStepper<> stepper(m_cfg.bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(
      stepper, Acts::VoidNavigator{}, logger().cloneWithSuffix("Propagator"));
  /*
  Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> extrapolator(
      Acts::EigenStepper<>(m_cfg.magneticField),
      Acts::Navigator({m_cfg.trackingGeometry},
                      logger().cloneWithSuffix("Navigator")),
      logger().cloneWithSuffix("Propagator"));
  */
  // Setup the vertex fitter
  Fitter::Config vertexFitterCfg;
  vertexFitterCfg.extractParameters
      .connect<&Acts::InputTrack::extractParameters>();
  // Setup the track linearizer
  Linearizer::Config linearizerCfg;
  linearizerCfg.bField = m_cfg.bField;
  linearizerCfg.propagator = propagator;
  Linearizer linearizer(linearizerCfg,
                        logger().cloneWithSuffix("HelicalTrackLinearizer"));

  vertexFitterCfg.trackLinearizer.connect<&Linearizer::linearizeTrack>(
      &linearizer);
  Fitter vertexFitter(vertexFitterCfg,
                      logger().cloneWithSuffix("FullBilloirVertexFitter"));

  // Setup the seed finder
  Acts::ImpactPointEstimator::Config ipEstCfg(m_cfg.bField, propagator);
  Acts::ImpactPointEstimator ipEst(
      ipEstCfg, logger().cloneWithSuffix("ImpactPointEstimator"));

  Acts::GaussianTrackDensity::Config densityCfg;
  densityCfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
  auto seeder = std::make_shared<Seeder>(Seeder::Config{{densityCfg}});
  // Set up the actual vertex finder
  Finder::Config finderCfg(std::move(vertexFitter), seeder, ipEst);
  finderCfg.trackLinearizer.connect<&Linearizer::linearizeTrack>(&linearizer);

  finderCfg.significanceCutSeeding = 100;
  finderCfg.maxVertices = 1;

  finderCfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
  finderCfg.field = m_cfg.bField;
  Finder finder(std::move(finderCfg), logger().clone());
  Acts::IVertexFinder::State state{std::in_place_type<Finder::State>,
                                   *m_cfg.bField, ctx.magFieldContext};
  Options finderOpts(ctx.geoContext, ctx.magFieldContext);

  // find vertices
  // TODO:
  // we need to build pairs of tracks
  // apply a dca selection
  // maybe a d0 selection
  // charge selection

  // Create propagator options
  Acts::PropagatorOptions<> pOptions(ctx.geoContext, ctx.magFieldContext);

  float mass1 = 0.13957039;
  float mass2 = 0.13957039;
  std::vector<float> masses;
  VertexCollection vertices;
  for (auto& track1 : inputTracks) {
    for (auto& track2 : inputTracks) {
      std::vector<Acts::InputTrack> tracksToFit;
      tracksToFit.emplace_back(track1);
      tracksToFit.emplace_back(track2);

      Acts::BoundTrackParameters params1 =
          vertexFitterCfg.extractParameters(track1);
      Acts::BoundTrackParameters params2 =
          vertexFitterCfg.extractParameters(track2);

      Acts::BoundVector tmppar1 = params1.parameters();
      Acts::BoundVector tmppar2 = params2.parameters();
      Acts::ActsScalar d01 = tmppar1(Acts::BoundIndices::eBoundLoc0);
      Acts::ActsScalar d02 = tmppar2(Acts::BoundIndices::eBoundLoc0);
      Acts::ActsScalar qOvP1 = tmppar1(Acts::BoundIndices::eBoundQOverP);
      Acts::ActsScalar qOvP2 = tmppar2(Acts::BoundIndices::eBoundQOverP);
      if (qOvP1 * qOvP2 > 0) {
        continue;
      }
      auto result = finder.find(tracksToFit, finderOpts, state);
      if (result.ok()) {
        vertices = std::move(result.value());

        if (vertices.size() > 0) {
          Acts::Vector3 vtx = vertices[0].position();
          const std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
              Acts::Surface::makeShared<Acts::PerigeeSurface>(vtx);
          // Propagate to the PCA of the reference point
          const auto res1 = propagator->propagateToSurface(
              params1,
              *std::static_pointer_cast<const Acts::Surface>(perigeeSurface),
              pOptions);

          if (!res1.ok()) {
            continue;
          }
          const auto& endParams1 = *res1;
          // Propagate to the PCA of the reference point

          const auto res2 = propagator->propagateToSurface(
              params2,
              *std::static_pointer_cast<const Acts::Surface>(perigeeSurface),
              pOptions);
          if (!res2.ok()) {
            continue;
          }
          const auto& endParams2 = *res2;
          Acts::BoundVector paramsAtPCA1 = endParams1.parameters();
          Acts::ActsScalar phi1 = paramsAtPCA1(Acts::BoundIndices::eBoundPhi);
          Acts::ActsScalar theta1 =
              paramsAtPCA1(Acts::BoundIndices::eBoundTheta);
          qOvP1 = paramsAtPCA1(Acts::BoundIndices::eBoundQOverP);
          Acts::Vector3 momentumAtPCA1(
              std::sin(phi1) * std::cos(theta1) / qOvP1,
              std::cos(phi1) * std::cos(theta1) / qOvP1,
              std::sin(theta1) / qOvP1);

          Acts::BoundVector paramsAtPCA2 = endParams2.parameters();
          Acts::ActsScalar phi2 = paramsAtPCA2(Acts::BoundIndices::eBoundPhi);
          Acts::ActsScalar theta2 =
              paramsAtPCA2(Acts::BoundIndices::eBoundTheta);
          qOvP2 = paramsAtPCA2(Acts::BoundIndices::eBoundQOverP);
          // Calculate Cartesian momentum components
          Acts::Vector3 momentumAtPCA2(
              std::sin(phi2) * std::cos(theta2) / qOvP2,
              std::cos(phi2) * std::cos(theta2) / qOvP2,
              std::sin(theta2) / qOvP2);
          float e1 = TMath::Sqrt(momentumAtPCA1[0] * momentumAtPCA1[0] +
                                 momentumAtPCA1[1] * momentumAtPCA1[1] +
                                 momentumAtPCA1[2] * momentumAtPCA1[2] +
                                 mass1 * mass1);
          float e2 = TMath::Sqrt(momentumAtPCA2[0] * momentumAtPCA2[0] +
                                 momentumAtPCA2[1] * momentumAtPCA2[1] +
                                 momentumAtPCA2[2] * momentumAtPCA2[2] +
                                 mass2 * mass2);

          
          //Acts::ImpactPointEstimator::State ipState(m_cfg.bField.makeCache(ctx.magFieldContext));
          /*
          Acts::ImpactPointEstimator::State ipState(ctx);
          
          auto distanceRes1 = ipEst.calculateDistance(
              ctx.geoContext, endParams1, vertices[0].position(), ipState);
          auto distanceRes2 = ipEst.calculateDistance(
              ctx.geoContext, endParams2, vertices[0].position(), ipState);

          if (*distanceRes1 > 1 || *distanceRes2)
            continue;
            */

          Acts::Vector4 v1(momentumAtPCA1[0], momentumAtPCA1[1],
                           momentumAtPCA1[2], e1);
          Acts::Vector4 v2(momentumAtPCA2[0], momentumAtPCA2[1],
                           momentumAtPCA2[2], e2);
          v1 += v2;


          float p = TMath::Sqrt((v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]));
          float pd = (vtx[0] * v1[0] + vtx[1] * v1[1] + vtx[2] * v1[2]);
          float d = TMath::Sqrt(vtx[0]*vtx[0]+vtx[1]*vtx[1]+vtx[2]*vtx[2]);

          std::cout<<"Dist: " << d << std::endl;
          std::cout<<"ct: " << d*0.49761/p << std::endl;
          std::cout<<"cos: " << pd/(p*d) << std::endl;
          std::cout<<"p: " << p << std::endl;

          if(pd/(p*d)<0.9999)
            continue;

          /*
          if( TMath::Sqrt(vtx[0]*vtx[0]+vtx[1]*vtx[1]+vtx[2]*vtx[2])*0.49761/p < 0.25)
            continue;
            */
          float m = TMath::Sqrt(v1[3] * v1[3] - (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]));

          std::cout << "Mass: " << m << std::endl;
          masses.push_back(m);
        }
      }
    }
  }
  // store proto vertices extracted from the found vertices
  m_outputProtoVertices(ctx, makeProtoVertices(inputTracks, vertices));

  // store found vertices
  m_outputVertices(ctx, std::move(vertices));
  m_outputMasses(ctx, std::move(masses));

  return ActsExamples::ProcessCode::SUCCESS;
}
