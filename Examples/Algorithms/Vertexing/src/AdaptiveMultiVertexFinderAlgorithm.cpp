// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
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
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include "VertexingHelpers.hpp"

ActsExamples::AdaptiveMultiVertexFinderAlgorithm::
    AdaptiveMultiVertexFinderAlgorithm(const Config& cfg,
                                       Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("AdaptiveMultiVertexFinder", lvl),
      m_cfg(cfg) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing input track parameters collection");
  }
  if (m_cfg.outputProtoVertices.empty()) {
    throw std::invalid_argument("Missing output proto vertices collection");
  }
  if (m_cfg.outputVertices.empty()) {
    throw std::invalid_argument("Missing output vertices collection");
  }
}

ActsExamples::ProcessCode
ActsExamples::AdaptiveMultiVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // retrieve input tracks and convert into the expected format
  const auto& inputTrackParameters =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  const auto& inputTrackPointers =
      makeTrackParametersPointerContainer(inputTrackParameters);

  //////////////////////////////////////////////
  /* Full tutorial example code for reference */
  //////////////////////////////////////////////

  // Set up EigenStepper
  Acts::EigenStepper<> stepper(m_cfg.bField);

  // Set up the propagator
  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
  auto propagator = std::make_shared<Propagator>(stepper);

  // Set up ImpactPointEstimator
  using IPEstimator =
      Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;
  IPEstimator::Config ipEstimatorCfg(m_cfg.bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  // Set up the helical track linearizer
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  Linearizer::Config ltConfig(m_cfg.bField, propagator);
  Linearizer linearizer(ltConfig);

  // Set up deterministic annealing with user-defined temperatures
  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  Acts::AnnealingUtility::Config annealingConfig(temperatures);
  Acts::AnnealingUtility annealingUtility(annealingConfig);

  // Set up the vertex fitter with user-defined annealing
  using Fitter =
      Acts::AdaptiveMultiVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  Fitter::Config fitterCfg(ipEstimator);
  fitterCfg.annealingTool = annealingUtility;
  Fitter fitter(fitterCfg);

  // Set up the vertex seed finder
  using SeedFinder = Acts::TrackDensityVertexFinder<
      Fitter, Acts::GaussianTrackDensity<Acts::BoundTrackParameters>>;
  SeedFinder seedFinder;

  // The vertex finder type
  using Finder = Acts::AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              linearizer, m_cfg.bField);
  // We do not want to use a beamspot constraint here
  finderConfig.useBeamSpotConstraint = false;

  // Instantiate the finder
  Finder finder(finderConfig);
  // The vertex finder state
  Finder::State state;

  // Default vertexing options, this is where e.g. a constraint could be set
  using VertexingOptions = Acts::VertexingOptions<Acts::BoundTrackParameters>;
  VertexingOptions finderOpts(ctx.geoContext, ctx.magFieldContext);

  // find vertices
  auto result = finder.find(inputTrackPointers, finderOpts, state);
  if (not result.ok()) {
    ACTS_ERROR("Error in vertex finder: " << result.error().message());
    return ProcessCode::ABORT;
  }
  auto vertices = *result;

  // show some debug output
  ACTS_INFO("Found " << vertices.size() << " vertices in event");
  for (const auto& vtx : vertices) {
    ACTS_DEBUG("Found vertex at " << vtx.fullPosition().transpose() << " with "
                                  << vtx.tracks().size() << " tracks.");
  }

  // store proto vertices extracted from the found vertices
  ctx.eventStore.add(m_cfg.outputProtoVertices,
                     makeProtoVertices(inputTrackParameters, vertices));

  std::vector<Acts::Vertex<Acts::BoundTrackParameters>> verticesOut = vertices;
  // store found vertices
  ctx.eventStore.add(m_cfg.outputVertices, std::move(vertices));

  return ActsExamples::ProcessCode::SUCCESS;
}
