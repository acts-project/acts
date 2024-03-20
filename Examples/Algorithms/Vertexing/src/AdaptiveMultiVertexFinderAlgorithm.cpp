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
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
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
      m_cfg(config),
      m_propagator{[&]() {
        // Set up EigenStepper
        Acts::EigenStepper<> stepper(m_cfg.bField);

        // Set up the propagator
        return std::make_shared<Propagator>(stepper);
      }()},
      m_ipEstimator{[&]() {
        // Set up ImpactPointEstimator
        Acts::ImpactPointEstimator::Config ipEstimatorCfg(m_cfg.bField,
                                                          m_propagator);
        return Acts::ImpactPointEstimator(
            ipEstimatorCfg, logger().cloneWithSuffix("ImpactPointEstimator"));
      }()},
      m_linearizer{[&] {
        // Set up the helical track linearizer
        Linearizer::Config ltConfig;
        ltConfig.bField = m_cfg.bField;
        ltConfig.propagator = m_propagator;
        return Linearizer(ltConfig,
                          logger().cloneWithSuffix("HelicalTrackLinearizer"));
      }()},
      m_vertexFinder{makeVertexFinder()} {
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

auto ActsExamples::AdaptiveMultiVertexFinderAlgorithm::makeVertexFinder() const
    -> Acts::AdaptiveMultiVertexFinder {
  std::shared_ptr<const Acts::IVertexFinder> seedFinder;
  if (m_cfg.seedFinder == SeedFinder::GaussianSeeder) {
    using Seeder = Acts::TrackDensityVertexFinder;
    Acts::GaussianTrackDensity::Config trkDensityCfg;
    trkDensityCfg.extractParameters
        .connect<&Acts::InputTrack::extractParameters>();
    seedFinder = std::make_shared<Seeder>(Seeder::Config{trkDensityCfg});
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
    using Seeder = Acts::AdaptiveGridDensityVertexFinder;
    Seeder::Config seederConfig(trkDensity);
    seederConfig.extractParameters
        .connect<&Acts::InputTrack::extractParameters>();
    seedFinder = std::make_shared<Seeder>(seederConfig);
  } else {
    throw std::invalid_argument("Unknown seed finder");
  }

  // Set up deterministic annealing with user-defined temperatures
  Acts::AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = {1.};
  Acts::AnnealingUtility annealingUtility(annealingConfig);

  // Set up the vertex fitter with user-defined annealing
  Fitter::Config fitterCfg(m_ipEstimator);
  fitterCfg.annealingTool = annealingUtility;
  fitterCfg.minWeight = 0.001;
  fitterCfg.doSmoothing = true;
  fitterCfg.useTime = m_cfg.useTime;
  fitterCfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
  fitterCfg.trackLinearizer.connect<&Linearizer::linearizeTrack>(&m_linearizer);
  Fitter fitter(std::move(fitterCfg),
                logger().cloneWithSuffix("AdaptiveMultiVertexFitter"));

  Acts::AdaptiveMultiVertexFinder::Config finderConfig(
      std::move(fitter), seedFinder, m_ipEstimator, m_cfg.bField);
  // Set the initial variance of the 4D vertex position. Since time is on a
  // numerical scale, we have to provide a greater value in the corresponding
  // dimension.
  finderConfig.initialVariances << 1e+2, 1e+2, 1e+2, 1e+8;
  finderConfig.tracksMaxZinterval = 1. * Acts::UnitConstants::mm;
  finderConfig.maxIterations = 200;
  finderConfig.useTime = m_cfg.useTime;
  // 5 corresponds to a p-value of ~0.92 using `chi2(x=5,ndf=2)`
  finderConfig.tracksMaxSignificance = 5;
  // This should be used consistently with and without time
  finderConfig.doFullSplitting = false;
  // 3 corresponds to a p-value of ~0.92 using `chi2(x=3,ndf=1)`
  finderConfig.maxMergeVertexSignificance = 3;
  if (m_cfg.useTime) {
    // When using time, we have an extra contribution to the chi2 by the time
    // coordinate. We thus need to increase tracksMaxSignificance (i.e., the
    // maximum chi2 that a track can have to be associated with a vertex).
    // Using the same p-value for 3 dof instead of 2.
    // 6.7 corresponds to a p-value of ~0.92 using `chi2(x=6.7,ndf=3)`
    finderConfig.tracksMaxSignificance = 6.7;
    // Using the same p-value for 2 dof instead of 1.
    // 5 corresponds to a p-value of ~0.92 using `chi2(x=5,ndf=2)`
    finderConfig.maxMergeVertexSignificance = 5;
  }
  finderConfig.extractParameters
      .template connect<&Acts::InputTrack::extractParameters>();

  // Instantiate the finder
  return Acts::AdaptiveMultiVertexFinder(std::move(finderConfig),
                                         logger().clone());
}

ActsExamples::ProcessCode
ActsExamples::AdaptiveMultiVertexFinderAlgorithm::execute(
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

  //////////////////////////////////////////////
  /* Full tutorial example code for reference */
  //////////////////////////////////////////////

  // The vertex finder state
  auto state = m_vertexFinder.makeState(ctx.magFieldContext);

  // Default vertexing options, this is where e.g. a constraint could be set
  Options finderOpts(ctx.geoContext, ctx.magFieldContext);

  VertexCollection vertices;

  if (inputTrackParameters.empty()) {
    ACTS_DEBUG("Empty track parameter collection found, skipping vertexing");
  } else {
    ACTS_DEBUG("Have " << inputTrackParameters.size()
                       << " input track parameters, running vertexing");
    // find vertices
    auto result = m_vertexFinder.find(inputTracks, finderOpts, state);

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
