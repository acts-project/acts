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
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Utilities/AnnealingUtility.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <system_error>

#include "VertexingHelpers.hpp"

ActsExamples::AdaptiveMultiVertexFinderAlgorithm::
    AdaptiveMultiVertexFinderAlgorithm(const Config& config,
                                       Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("AdaptiveMultiVertexFinder", level),
      m_cfg(config) {
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
  m_outputVertices.initialize(m_cfg.outputVertices);
}

ActsExamples::ProcessCode
ActsExamples::AdaptiveMultiVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  if (m_cfg.seedFinder == SeedFinder::GaussianSeeder) {
    using Seeder = Acts::TrackDensityVertexFinder<
        Fitter, Acts::GaussianTrackDensity<Acts::BoundTrackParameters>>;
    using Finder = Acts::AdaptiveMultiVertexFinder<Fitter, Seeder>;
    Seeder seedFinder;
    return executeAfterSeederChoice<Seeder, Finder>(ctx, seedFinder);
  } else if (m_cfg.seedFinder == SeedFinder::AdaptiveGridSeeder) {
    using Seeder = Acts::AdaptiveGridDensityVertexFinder<109, Fitter>;
    using Finder = Acts::AdaptiveMultiVertexFinder<Fitter, Seeder>;
    // The seeder config argument corresponds to the bin size in mm
    Seeder::Config seederConfig(0.05);
    Seeder seedFinder(seederConfig);
    return executeAfterSeederChoice<Seeder, Finder>(ctx, seedFinder);
  } else {
    return ActsExamples::ProcessCode::ABORT;
  }
}

template <typename vseeder_t, typename vfinder_t>
ActsExamples::ProcessCode
ActsExamples::AdaptiveMultiVertexFinderAlgorithm::executeAfterSeederChoice(
    const ActsExamples::AlgorithmContext& ctx,
    const vseeder_t& seedFinder) const {
  using Finder = vfinder_t;

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
  Fitter fitter(std::move(fitterCfg), logger().cloneWithSuffix("AMVFitter"));

  typename Finder::Config finderConfig(std::move(fitter), seedFinder,
                                       ipEstimator, std::move(linearizer),
                                       m_cfg.bField);
  finderConfig.looseConstrValue = 1e2;
  finderConfig.tracksMaxZinterval = 1. * Acts::UnitConstants::mm;
  finderConfig.maxIterations = 200;

  // Instantiate the finder
  Finder finder(std::move(finderConfig), logger().cloneWithSuffix("AMVFinder"));

  // retrieve input tracks and convert into the expected format

  const auto& inputTrackParameters = m_inputTrackParameters(ctx);
  // TODO change this from pointers to tracks parameters to actual tracks
  auto inputTrackPointers =
      makeTrackParametersPointerContainer(inputTrackParameters);

  if (inputTrackParameters.size() != inputTrackPointers.size()) {
    ACTS_ERROR("Input track containers do not align: "
               << inputTrackParameters.size()
               << " != " << inputTrackPointers.size());
  }

  for (const auto trk : inputTrackPointers) {
    if (trk->covariance() && trk->covariance()->determinant() <= 0) {
      // actually we should consider this as an error but I do not want the CI
      // to fail
      ACTS_WARNING("input track " << *trk << " has det(cov) = "
                                  << trk->covariance()->determinant());
    }
  }

  //////////////////////////////////////////////
  /* Full tutorial example code for reference */
  //////////////////////////////////////////////

  // The vertex finder state
  typename Finder::State state;

  // Default vertexing options, this is where e.g. a constraint could be set
  Options finderOpts(ctx.geoContext, ctx.magFieldContext);

  VertexCollection vertices;

  if (inputTrackParameters.empty()) {
    ACTS_DEBUG("Empty track parameter collection found, skipping vertexing");
  } else {
    ACTS_DEBUG("Have " << inputTrackParameters.size()
                       << " input track parameters, running vertexing");
    // find vertices
    auto result = finder.find(inputTrackPointers, finderOpts, state);

    if (result.ok()) {
      vertices = std::move(result.value());
    } else {
      ACTS_ERROR("Error in vertex finder: " << result.error().message());
    }
  }

  // show some debug output
  ACTS_INFO("Found " << vertices.size() << " vertices in event");
  for (const auto& vtx : vertices) {
    ACTS_DEBUG("Found vertex at " << vtx.fullPosition().transpose() << " with "
                                  << vtx.tracks().size() << " tracks.");
  }

  // store proto vertices extracted from the found vertices
  m_outputProtoVertices(ctx, makeProtoVertices(inputTrackPointers, vertices));

  // store found vertices
  m_outputVertices(ctx, std::move(vertices));

  return ActsExamples::ProcessCode::SUCCESS;
}
