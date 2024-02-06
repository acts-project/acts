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
#include "Acts/Vertexing/AdaptiveGridTrackDensity.hpp"
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
    // Set up track density used during vertex seeding
    Acts::AdaptiveGridTrackDensity::Config trkDensityCfg;
    // Bin extent in z-direction
    trkDensityCfg.spatialBinExtent = 15. * Acts::UnitConstants::um;
    // Bin extent in t-direction
    trkDensityCfg.temporalBinExtent = 19. * Acts::UnitConstants::mm;
    trkDensityCfg.useTime = m_cfg.useTime;
    Acts::AdaptiveGridTrackDensity trkDensity(trkDensityCfg);

    // Set up vertex seeder and finder
    using Seeder = Acts::AdaptiveGridDensityVertexFinder<Fitter>;
    using Finder = Acts::AdaptiveMultiVertexFinder<Fitter, Seeder>;
    Seeder::Config seederConfig(trkDensity);
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
  IPEstimator ipEstimator(ipEstimatorCfg,
                          logger().cloneWithSuffix("ImpactPointEstimator"));

  // Set up the helical track linearizer
  Linearizer::Config ltConfig(m_cfg.bField, propagator);
  Linearizer linearizer(ltConfig,
                        logger().cloneWithSuffix("HelicalTrackLinearizer"));

  // Set up deterministic annealing with user-defined temperatures
  Acts::AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = {1.};
  Acts::AnnealingUtility annealingUtility(annealingConfig);

  // Set up the vertex fitter with user-defined annealing
  Fitter::Config fitterCfg(ipEstimator);
  fitterCfg.annealingTool = annealingUtility;
  fitterCfg.minWeight = 0.001;
  fitterCfg.doSmoothing = true;
  fitterCfg.useTime = m_cfg.useTime;
  Fitter fitter(std::move(fitterCfg),
                logger().cloneWithSuffix("AdaptiveMultiVertexFitter"));

  typename Finder::Config finderConfig(std::move(fitter), std::move(seedFinder),
                                       ipEstimator, std::move(linearizer),
                                       m_cfg.bField);
  // Set the initial variance of the 4D vertex position. Since time is on a
  // numerical scale, we have to provide a greater value in the corresponding
  // dimension.
  finderConfig.initialVariances << 1e+2, 1e+2, 1e+2, 1e+8;
  finderConfig.tracksMaxZinterval = 1. * Acts::UnitConstants::mm;
  finderConfig.maxIterations = 200;
  finderConfig.useTime = m_cfg.useTime;
  if (m_cfg.useTime) {
    // When using time, we have an extra contribution to the chi2 by the time
    // coordinate. We thus need to increase tracksMaxSignificance (i.e., the
    // maximum chi2 that a track can have to be associated with a vertex).
    finderConfig.tracksMaxSignificance = 7.5;
    // Check if vertices are merged in space and time
    // TODO rename do3dSplitting -> doFullSplitting
    finderConfig.do3dSplitting = true;
    // Reset the maximum significance that two vertices can have before they are
    // considered as merged. The default value 3 is tuned for comparing the
    // vertices' z-coordinates. Since we consider 4 dimensions here, we need to
    // multiply the value by 4 and thus we set it to 3 * 4 = 12.
    finderConfig.maxMergeVertexSignificance = 12.;
  }

  // Instantiate the finder
  Finder finder(std::move(finderConfig), logger().clone());

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
    auto result = finder.find(inputTracks, finderOpts, state);

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
  m_outputProtoVertices(ctx, makeProtoVertices(inputTracks, vertices));

  // store found vertices
  m_outputVertices(ctx, std::move(vertices));

  return ActsExamples::ProcessCode::SUCCESS;
}
