// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Vertexing/VertexFitAlgorithm.hpp"

#include <iostream>

#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/TruthTracking/VertexAndTracks.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

FWE::VertexFitAlgorithm::VertexFitAlgorithm(const Config& cfg,
                                            Acts::Logging::Level level)
    : FW::BareAlgorithm("VertexFit", level), m_cfg(cfg) {}

/// @brief Algorithm that receives a set of tracks belonging to a common
/// vertex and fits the associated vertex to it
FW::ProcessCode FWE::VertexFitAlgorithm::execute(
    const FW::AlgorithmContext& ctx) const {
  using MagneticField = Acts::ConstantBField;
  using Stepper = Acts::EigenStepper<MagneticField>;
  using Propagator = Acts::Propagator<Stepper>;
  using PropagatorOptions = Acts::PropagatorOptions<>;
  using TrackParameters = Acts::BoundParameters;
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  using VertexFitter =
      Acts::FullBilloirVertexFitter<TrackParameters, Linearizer>;
  using VertexFitterOptions = Acts::VertexingOptions<TrackParameters>;

  // Setup the magnetic field
  MagneticField bField(m_cfg.bField);
  // Setup the propagator with void navigator
  auto propagator = std::make_shared<Propagator>(Stepper(bField));
  PropagatorOptions propagatorOpts(ctx.geoContext, ctx.magFieldContext);
  // Setup the vertex fitter
  VertexFitter::Config vertexFitterCfg;
  VertexFitter vertexFitter(vertexFitterCfg);
  VertexFitter::State state(ctx.magFieldContext);
  // Setup the linearizer
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  const auto& input = ctx.eventStore.get<std::vector<FW::VertexAndTracks>>(
      m_cfg.trackCollection);
  for (auto& vertexAndTracks : input) {
    const auto& inputTrackCollection = vertexAndTracks.tracks;

    // Only fit vertices where true vertex is indeed close to beam line
    // This removed secondaries and odd vertices that are not wanted in
    // this example
    if (std::sqrt(vertexAndTracks.vertex.position().x() *
                      vertexAndTracks.vertex.position().x() +
                  vertexAndTracks.vertex.position().y() *
                      vertexAndTracks.vertex.position().y()) >
        m_cfg.maxTransFitDistance) {
      continue;
    }

    std::vector<const Acts::BoundParameters*> inputTrackPtrCollection;
    for (const auto& trk : inputTrackCollection) {
      inputTrackPtrCollection.push_back(&trk);
    }

    Acts::Vertex<TrackParameters> fittedVertex;
    if (!m_cfg.doConstrainedFit) {
      if (inputTrackCollection.size() < 2) {
        continue;
      }
      // Vertex fitter options
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
      Acts::Vertex<TrackParameters> theConstraint;

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

    ACTS_INFO("Fitted Vertex: "
              << "(" << fittedVertex.position().x() << ","
              << fittedVertex.position().y() << ","
              << fittedVertex.position().z() << ")");
    ACTS_INFO("Truth Vertex: "
              << "(" << vertexAndTracks.vertex.position().x() << ","
              << vertexAndTracks.vertex.position().y() << ","
              << vertexAndTracks.vertex.position().z() << ")");
  }

  return FW::ProcessCode::SUCCESS;
}
