// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/VertexFitterAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldQuery.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFitInput.hpp"
#include "ActsExamples/EventData/Vertex.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <optional>
#include <ostream>
#include <stdexcept>

namespace ActsExamples {

VertexFitterAlgorithm::VertexFitterAlgorithm(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("VertexFit", std::move(logger)), m_cfg(cfg) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing input track parameter collection");
  }
  if (m_cfg.inputProtoVertices.empty()) {
    throw std::invalid_argument("Missing input proto vertices collection");
  }

  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_inputProtoVertices.initialize(m_cfg.inputProtoVertices);
  m_outputVertices.initialize(m_cfg.outputVertices);
}

ProcessCode VertexFitterAlgorithm::execute(const AlgorithmContext& ctx) const {
  using Propagator = Acts::Propagator<Acts::SympyStepper>;
  using PropagatorOptions = Propagator::Options<>;
  using Linearizer = Acts::HelicalTrackLinearizer;
  using VertexFitter = Acts::FullBilloirVertexFitter;

  // Set up SympyStepper
  Acts::SympyStepper stepper(m_cfg.bField);

  // Setup the propagator with void navigator
  auto propagator = std::make_shared<Propagator>(
      stepper, Acts::VoidNavigator{}, logger().cloneWithSuffix("Prop"));

  // Setup the linearizer
  Linearizer::Config ltConfig;
  ltConfig.propagator = propagator;
  Linearizer linearizer(ltConfig, logger().cloneWithSuffix("HelLin"));

  PropagatorOptions propagatorOpts(ctx.recoGeoContext, ctx.magFieldContext);
  Acts::MagneticFieldQuery fieldQuery(m_cfg.bField, ctx.magFieldContext);

  // Setup the vertex fitter
  VertexFitter::Config vertexFitterCfg;
  vertexFitterCfg.extractParameters
      .connect<&Acts::InputTrack::extractParameters>();
  vertexFitterCfg.fieldTrackLinearizer
      .connect<&Linearizer::linearizeTrackWithField>(&linearizer);
  VertexFitter vertexFitter(vertexFitterCfg);

  ACTS_VERBOSE("Read from '" << m_cfg.inputTrackParameters << "'");
  ACTS_VERBOSE("Read from '" << m_cfg.inputProtoVertices << "'");

  const auto& inputTrackParameters = m_inputTrackParameters(ctx);
  ACTS_VERBOSE("Have " << inputTrackParameters.size() << " track parameters");
  const auto& protoVertices = m_inputProtoVertices(ctx);
  ACTS_VERBOSE("Have " << protoVertices.size() << " proto vertices");

  std::vector<Acts::InputTrack> inputTracks;

  VertexContainer fittedVertices;

  for (const auto& protoVertex : protoVertices) {
    // un-constrained fit requires at least two tracks
    if ((!m_cfg.doConstrainedFit) && (protoVertex.size() < 2)) {
      ACTS_INFO(
          "Skip un-constrained vertex fit on proto-vertex with less than two "
          "tracks");
      continue;
    }

    // select input tracks for the input proto vertex
    inputTracks.clear();
    inputTracks.reserve(protoVertex.size());
    for (const auto& trackIdx : protoVertex) {
      if (trackIdx >= inputTrackParameters.size()) {
        ACTS_ERROR("track parameters " << trackIdx << " does not exist");
        continue;
      }

      inputTracks.emplace_back(&inputTrackParameters[trackIdx]);
    }

    if (!m_cfg.doConstrainedFit) {
      Acts::VertexFitInput input{inputTracks, Acts::Vector4::Zero(),
                                 std::nullopt};
      auto fitRes = vertexFitter.fit(input, ctx.recoGeoContext,
                                     ctx.magFieldContext, fieldQuery);
      if (fitRes.ok()) {
        fittedVertices.push_back(*fitRes);
      } else {
        ACTS_ERROR("Error in vertex fitter: " << fitRes.error().message());
      }
    } else {
      // Vertex constraint
      Acts::Vertex theConstraint;

      theConstraint.setFullCovariance(m_cfg.constraintCov);
      theConstraint.setFullPosition(m_cfg.constraintPos);

      Acts::VertexFitInput input{inputTracks, theConstraint.fullPosition(),
                                 theConstraint};
      auto fitRes = vertexFitter.fit(input, ctx.recoGeoContext,
                                     ctx.magFieldContext, fieldQuery);
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

  m_outputVertices(ctx, std::move(fittedVertices));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
