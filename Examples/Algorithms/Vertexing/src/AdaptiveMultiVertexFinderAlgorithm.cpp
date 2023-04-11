// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <chrono>
#include <memory>

#include "VertexingHelpers.hpp"

ActsExamples::AdaptiveMultiVertexFinderAlgorithm::
    AdaptiveMultiVertexFinderAlgorithm(const Config& config,
                                       Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("AdaptiveMultiVertexFinder", level),
      m_cfg(config) {
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
  if (m_cfg.outputTime.empty()) {
    throw std::invalid_argument("Missing output reconstruction time");
  }

  m_inputTrackParameters.maybeInitialize(m_cfg.inputTrackParameters);
  m_inputTrajectories.maybeInitialize(m_cfg.inputTrajectories);

  m_outputProtoVertices.initialize(m_cfg.outputProtoVertices);
  m_outputVertices.initialize(m_cfg.outputVertices);
  m_outputTime.initialize(m_cfg.outputTime);
}

ActsExamples::ProcessCode
ActsExamples::AdaptiveMultiVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Set up EigenStepper
  Acts::EigenStepper<> stepper(m_cfg.bField);

  // Set up the propagator
  auto propagator = std::make_shared<Propagator>(stepper);

  // Set up ImpactPointEstimator
  IPEstimator::Config ipEstimatorCfg(m_cfg.bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  // Set up the helical track linearizer
  Linearizer::Config ltConfig(m_cfg.bField, propagator);
  Linearizer linearizer(ltConfig, logger().cloneWithSuffix("HelLin"));

  // Set up deterministic annealing with user-defined temperatures
  Acts::AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = {1.};
  Acts::AnnealingUtility annealingUtility(annealingConfig);

  // Set up the vertex fitter with user-defined annealing
  Fitter::Config fitterCfg(ipEstimator);
  fitterCfg.annealingTool = annealingUtility;
  fitterCfg.minWeight = 0.001;
  fitterCfg.doSmoothing = true;
  Fitter fitter(fitterCfg, logger().cloneWithSuffix("AMVFitter"));

  // Set up the vertex seed finder
  SeedFinder seedFinder;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              std::move(linearizer), m_cfg.bField);
  // We do not want to use a beamspot constraint here
  finderConfig.useBeamSpotConstraint = false;
  finderConfig.tracksMaxZinterval = 1. * Acts::UnitConstants::mm;

  // Instantiate the finder
  Finder finder(finderConfig, logger().cloneWithSuffix("AMVFinder"));

  // retrieve input tracks and convert into the expected format

  auto [inputTrackParameters, inputTrackPointers] =
      makeParameterContainers(ctx, m_inputTrackParameters, m_inputTrajectories);

  //////////////////////////////////////////////
  /* Full tutorial example code for reference */
  //////////////////////////////////////////////

  // The vertex finder state
  Finder::State state;

  // Default vertexing options, this is where e.g. a constraint could be set
  using VertexingOptions = Acts::VertexingOptions<Acts::BoundTrackParameters>;
  VertexingOptions finderOpts(ctx.geoContext, ctx.magFieldContext);

  VertexCollection vertices;

  auto t1 = std::chrono::high_resolution_clock::now();

  if (inputTrackParameters.empty()) {
    ACTS_DEBUG("Empty track parameter collection found, skipping vertexing");
  } else {
    ACTS_DEBUG("Have " << inputTrackParameters.size()
                       << " input track parameters, running vertexing");
    // find vertices and measure elapsed time
    auto result = finder.find(inputTrackPointers, finderOpts, state);

    if (result.ok()) {
      vertices = std::move(result.value());
    } else {
      ACTS_ERROR("Error in vertex finder: " << result.error().message());
    }
  }

  auto t2 = std::chrono::high_resolution_clock::now();

  // show some debug output
  ACTS_INFO("Found " << vertices.size() << " vertices in event");
  for (const auto& vtx : vertices) {
    ACTS_DEBUG("Found vertex at " << vtx.fullPosition().transpose() << " with "
                                  << vtx.tracks().size() << " tracks.");
  }

  // store proto vertices extracted from the found vertices
  m_outputProtoVertices(ctx, makeProtoVertices(inputTrackParameters, vertices));

  // store found vertices
  m_outputVertices(ctx, std::move(vertices));

  // time in milliseconds
  int timeMS =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  // store reconstruction time
  m_outputTime(ctx,
               std::move(timeMS));  // NOLINT(performance-move-const-arg)

  return ActsExamples::ProcessCode::SUCCESS;
}
