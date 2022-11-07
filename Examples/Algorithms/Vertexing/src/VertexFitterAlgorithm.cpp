// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/VertexFitterAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

ActsExamples::VertexFitterAlgorithm::VertexFitterAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("VertexFit", lvl), m_cfg(cfg) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing input track parameters collection");
  }
  if (m_cfg.inputProtoVertices.empty()) {
    throw std::invalid_argument("Missing input proto vertices collection");
  }
}

ActsExamples::ProcessCode ActsExamples::VertexFitterAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
  using PropagatorOptions = Acts::PropagatorOptions<>;
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  using VertexFitter =
      Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  using VertexFitterOptions =
      Acts::VertexingOptions<Acts::BoundTrackParameters>;

  // Set up EigenStepper
  Acts::EigenStepper<> stepper(m_cfg.bField);

  // Setup the propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);
  PropagatorOptions propagatorOpts(ctx.geoContext, ctx.magFieldContext,
                                   Acts::LoggerWrapper{logger()});
  // Setup the vertex fitter
  VertexFitter::Config vertexFitterCfg;
  VertexFitter vertexFitter(vertexFitterCfg);
  VertexFitter::State state(m_cfg.bField->makeCache(ctx.magFieldContext));
  // Setup the linearizer
  Linearizer::Config ltConfig(m_cfg.bField, propagator);
  Linearizer linearizer(ltConfig);

  ACTS_VERBOSE("Read from '" << m_cfg.inputTrackParameters << "'");
  ACTS_VERBOSE("Read from '" << m_cfg.inputProtoVertices << "'");

  const auto& trackParameters =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  ACTS_VERBOSE("Have " << trackParameters.size() << " track parameters");
  const auto& protoVertices =
      ctx.eventStore.get<ProtoVertexContainer>(m_cfg.inputProtoVertices);
  ACTS_VERBOSE("Have " << protoVertices.size() << " proto vertices");

  std::vector<const Acts::BoundTrackParameters*> inputTrackPtrCollection;

  std::vector<Acts::Vertex<Acts::BoundTrackParameters>> fittedVertices;

  for (const auto& protoVertex : protoVertices) {
    // un-constrained fit requires at least two tracks
    if ((not m_cfg.doConstrainedFit) and (protoVertex.size() < 2)) {
      ACTS_INFO(
          "Skip un-constrained vertex fit on proto-vertex with less than two "
          "tracks");
      continue;
    }

    // select input tracks for the input proto vertex
    inputTrackPtrCollection.clear();
    inputTrackPtrCollection.reserve(protoVertex.size());
    for (const auto& trackIdx : protoVertex) {
      if (trackIdx >= trackParameters.size()) {
        ACTS_ERROR("track parameters " << trackIdx << " does not exist");
        continue;
      }

      inputTrackPtrCollection.push_back(&trackParameters[trackIdx]);
    }

    if (!m_cfg.doConstrainedFit) {
      VertexFitterOptions vfOptions(ctx.geoContext, ctx.magFieldContext);

      auto fitRes = vertexFitter.fit(inputTrackPtrCollection, linearizer,
                                     vfOptions, state);
      if (fitRes.ok()) {
        fittedVertices.push_back(*fitRes);
      } else {
        ACTS_ERROR("Error in vertex fitter: " << fitRes.error().message());
      }
    } else {
      // Vertex constraint
      Acts::Vertex<Acts::BoundTrackParameters> theConstraint;

      theConstraint.setFullCovariance(m_cfg.constraintCov);
      theConstraint.setFullPosition(m_cfg.constraintPos);

      // Vertex fitter options
      VertexFitterOptions vfOptionsConstr(ctx.geoContext, ctx.magFieldContext,
                                          theConstraint);

      auto fitRes = vertexFitter.fit(inputTrackPtrCollection, linearizer,
                                     vfOptionsConstr, state);
      if (fitRes.ok()) {
        fittedVertices.push_back(*fitRes);
      } else {
        ACTS_ERROR(
            "Error in constrained vertex fitter: " << fitRes.error().message());
      }
    }

    if (fittedVertices.empty()) {
      ACTS_DEBUG("No fitted vertex");
    } else {
      ACTS_DEBUG("Fitted Vertex "
                 << fittedVertices.back().fullPosition().transpose());
      ACTS_DEBUG(
          "Tracks at fitted Vertex: " << fittedVertices.back().tracks().size());
    }
  }

  ctx.eventStore.add(m_cfg.outputVertices, std::move(fittedVertices));
  return ProcessCode::SUCCESS;
}
