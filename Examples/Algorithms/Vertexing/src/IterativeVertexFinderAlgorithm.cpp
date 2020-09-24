// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "Acts/Vertexing/ZScanVertexFinder.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include "VertexingHelpers.hpp"

ActsExamples::IterativeVertexFinderAlgorithm::IterativeVertexFinderAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("IterativeVertexFinder", lvl), m_cfg(cfg) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing input track parameters collection");
  }
  if (m_cfg.outputProtoVertices.empty()) {
    throw std::invalid_argument("Missing output proto vertices collection");
  }
}

ActsExamples::ProcessCode ActsExamples::IterativeVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // retrieve input tracks and convert into the expected format
  const auto& inputTrackParameters =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  const auto& inputTrackPointers =
      makeTrackParametersPointerContainer(inputTrackParameters);

  using MagneticField = Acts::ConstantBField;
  using Stepper = Acts::EigenStepper<MagneticField>;
  using Propagator = Acts::Propagator<Stepper>;
  using PropagatorOptions = Acts::PropagatorOptions<>;
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  using VertexFitter =
      Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  using ImpactPointEstimator =
      Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;
  using VertexSeeder = Acts::ZScanVertexFinder<VertexFitter>;
  using VertexFinder = Acts::IterativeVertexFinder<VertexFitter, VertexSeeder>;
  using VertexFinderOptions =
      Acts::VertexingOptions<Acts::BoundTrackParameters>;

  // Set up the magnetic field
  MagneticField bField(m_cfg.bField);
  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(Stepper(bField));
  PropagatorOptions propagatorOpts(ctx.geoContext, ctx.magFieldContext,
                                   Acts::LoggerWrapper{logger()});
  // Setup the vertex fitter
  VertexFitter::Config vertexFitterCfg;
  VertexFitter vertexFitter(std::move(vertexFitterCfg));
  // Setup the track linearizer
  Linearizer::Config linearizerCfg(bField, propagator);
  Linearizer linearizer(std::move(linearizerCfg));
  // Setup the seed finder
  ImpactPointEstimator::Config ipEstCfg(bField, propagator);
  ImpactPointEstimator ipEst(std::move(ipEstCfg));
  VertexSeeder::Config seederCfg(ipEst);
  VertexSeeder seeder(std::move(seederCfg));
  // Set up the actual vertex finder
  VertexFinder::Config finderCfg(std::move(vertexFitter), std::move(linearizer),
                                 std::move(seeder), ipEst);
  finderCfg.maxVertices = 200;
  finderCfg.reassignTracksAfterFirstFit = true;
  VertexFinder finder(finderCfg);
  VertexFinder::State state(ctx.magFieldContext);
  VertexFinderOptions finderOpts(ctx.geoContext, ctx.magFieldContext);

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
    ACTS_INFO("Found vertex at " << vtx.fullPosition().transpose() << " with "
                                 << vtx.tracks().size() << " tracks.");
  }

  // store proto vertices extracted from the found vertices
  ctx.eventStore.add(m_cfg.outputProtoVertices,
                     makeProtoVertices(inputTrackParameters, vertices));

  return ActsExamples::ProcessCode::SUCCESS;
}
