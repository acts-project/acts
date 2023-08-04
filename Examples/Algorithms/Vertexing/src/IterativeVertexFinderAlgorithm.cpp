// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/detail/VoidPropagatorComponents.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <chrono>
#include <ostream>
#include <stdexcept>
#include <system_error>

#include "VertexingHelpers.hpp"

ActsExamples::IterativeVertexFinderAlgorithm::IterativeVertexFinderAlgorithm(
    const Config& config, Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("IterativeVertexFinder", level), m_cfg(config) {
  if (m_cfg.inputTrackParameters.empty() == m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument(
        "You have to either provide track parameters or trajectories");
  }
  if (m_cfg.outputProtoVertices.empty()) {
    throw std::invalid_argument("Missing output proto vertices collection");
  }
  if (m_cfg.outputVertices.empty()) {
    throw std::invalid_argument("Missing output vertices collection");
  }

  m_inputTrackParameters.maybeInitialize(m_cfg.inputTrackParameters);
  m_inputTrajectories.maybeInitialize(m_cfg.inputTrajectories);

  m_outputProtoVertices.initialize(m_cfg.outputProtoVertices);
  m_outputVertices.initialize(m_cfg.outputVertices);
}

ActsExamples::ProcessCode ActsExamples::IterativeVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // retrieve input tracks and convert into the expected format

  auto [inputTrackParameters, inputTrackPointers] =
      makeParameterContainers(ctx, m_inputTrackParameters, m_inputTrajectories);

  if (inputTrackParameters.size() != inputTrackPointers.size()) {
    ACTS_ERROR("Input track containers do not align: "
               << inputTrackParameters.size()
               << " != " << inputTrackPointers.size());
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
      stepper, Acts::detail::VoidNavigator{}, logger().cloneWithSuffix("Prop"));
  // Setup the vertex fitter
  Fitter::Config vertexFitterCfg;
  Fitter vertexFitter(vertexFitterCfg);
  // Setup the track linearizer
  Linearizer::Config linearizerCfg(m_cfg.bField, propagator);
  Linearizer linearizer(linearizerCfg, logger().cloneWithSuffix("HelLin"));
  // Setup the seed finder
  IPEstimator::Config ipEstCfg(m_cfg.bField, propagator);
  IPEstimator ipEst(ipEstCfg);
  Seeder seeder;
  // Set up the actual vertex finder
  Finder::Config finderCfg(vertexFitter, std::move(linearizer),
                           std::move(seeder), ipEst);
  finderCfg.maxVertices = 200;
  finderCfg.reassignTracksAfterFirstFit = false;
  Finder finder(finderCfg, logger().cloneWithSuffix("Finder"));
  Finder::State state(*m_cfg.bField, ctx.magFieldContext);
  Options finderOpts(ctx.geoContext, ctx.magFieldContext);

  // find vertices
  auto result = finder.find(inputTrackPointers, finderOpts, state);

  VertexCollection vertices;
  if (result.ok()) {
    vertices = std::move(result.value());
  } else {
    ACTS_ERROR("Error in vertex finder: " << result.error().message());
  }

  // show some debug output
  ACTS_INFO("Found " << vertices.size() << " vertices in event");
  for (const auto& vtx : vertices) {
    ACTS_INFO("Found vertex at " << vtx.fullPosition().transpose() << " with "
                                 << vtx.tracks().size() << " tracks.");
  }

  // store proto vertices extracted from the found vertices
  m_outputProtoVertices(ctx, makeProtoVertices(inputTrackParameters, vertices));

  // store found vertices
  m_outputVertices(ctx, std::move(vertices));

  return ActsExamples::ProcessCode::SUCCESS;
}
