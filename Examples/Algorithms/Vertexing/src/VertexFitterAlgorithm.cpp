// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/VertexFitterAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
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
  using MagneticField = Acts::ConstantBField;
  using Stepper = Acts::EigenStepper<MagneticField>;
  using Propagator = Acts::Propagator<Stepper>;
  using PropagatorOptions = Acts::PropagatorOptions<>;
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  using VertexFitter =
      Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  using VertexFitterOptions =
      Acts::VertexingOptions<Acts::BoundTrackParameters>;

  // Setup the magnetic field
  MagneticField bField(m_cfg.bField);
  // Setup the propagator with void navigator
  auto propagator = std::make_shared<Propagator>(Stepper(bField));
  PropagatorOptions propagatorOpts(ctx.geoContext, ctx.magFieldContext,
                                   Acts::LoggerWrapper{logger()});
  // Setup the vertex fitter
  VertexFitter::Config vertexFitterCfg;
  VertexFitter vertexFitter(vertexFitterCfg);
  VertexFitter::State state(ctx.magFieldContext);
  // Setup the linearizer
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  const auto& trackParameters =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  const auto& protoVertices =
      ctx.eventStore.get<ProtoVertexContainer>(m_cfg.inputProtoVertices);
  std::vector<const Acts::BoundTrackParameters*> inputTrackPtrCollection;

  for (const auto& protoVertex : protoVertices) {
    // un-constrained fit requires at least two tracks
    if ((not m_cfg.doConstrainedFit) and (protoVertex.size() < 2)) {
      ACTS_WARNING(
          "Skip un-constrained vertex fit on proto-vertex with less than two "
          "tracks");
      continue;
    }

    // select input tracks for the input proto vertex
    inputTrackPtrCollection.clear();
    inputTrackPtrCollection.reserve(protoVertex.size());
    for (const auto& trackIdx : protoVertex) {
      inputTrackPtrCollection.push_back(&trackParameters[trackIdx]);
    }

    Acts::Vertex<Acts::BoundTrackParameters> fittedVertex;
    if (!m_cfg.doConstrainedFit) {
      VertexFitterOptions vfOptions(ctx.geoContext, ctx.magFieldContext);

      auto fitRes = vertexFitter.fit(inputTrackPtrCollection, linearizer,
                                     vfOptions, state);
      if (fitRes.ok()) {
        fittedVertex = *fitRes;
      } else {
        ACTS_ERROR("Error in vertex fit.");
        ACTS_ERROR(fitRes.error().message());
      }
    } else {
      // Vertex constraint
      Acts::Vertex<Acts::BoundTrackParameters> theConstraint;

      theConstraint.setCovariance(m_cfg.constraintCov);
      theConstraint.setPosition(m_cfg.constraintPos);

      // Vertex fitter options
      VertexFitterOptions vfOptionsConstr(ctx.geoContext, ctx.magFieldContext,
                                          theConstraint);

      auto fitRes = vertexFitter.fit(inputTrackPtrCollection, linearizer,
                                     vfOptionsConstr, state);
      if (fitRes.ok()) {
        fittedVertex = *fitRes;
      } else {
        ACTS_ERROR("Error in vertex fit with constraint.");
        ACTS_ERROR(fitRes.error().message());
      }
    }

    ACTS_INFO("Fitted Vertex " << fittedVertex.fullPosition().transpose());
  }
  return ProcessCode::SUCCESS;
}
