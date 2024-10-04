// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TransportParticles/TransportParticlesAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

#include "TRandom.h"

ActsExamples::TransportParticles::TransportParticles(
    ActsExamples::TransportParticles::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("TransportParticles", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackParameters) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackParameters.emplace_back(
        std::make_unique<ReadDataHandle<TrackParametersContainer>>(
            this, "inputTrackParameters#" +
                      std::to_string(m_inputTrackParameters.size())));
    handle->initialize(spName);
  }

  if (m_cfg.inputTrackContainer.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackContainer) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackContainer.emplace_back(
        std::make_unique<ReadDataHandle<ConstTrackContainer>>(
            this, "inputTrackContainer#" +
                      std::to_string(m_inputTrackContainer.size())));
    handle->initialize(spName);
  }

  m_inputVertices.initialize(m_cfg.inputVertices);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode ActsExamples::TransportParticles::execute(
    const AlgorithmContext& ctx) const {
  const auto& inputVertices = m_inputVertices(ctx);
  Acts::Vector3 propPoint = Acts::Vector3{m_cfg.px, m_cfg.py, m_cfg.pz};
  if (m_cfg.useRecVtx)
    propPoint = inputVertices[0].position();

  std::cout << "propPoint: " << propPoint[0] << " " << propPoint[1] << " "
            << propPoint[2] << std::endl;
  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(propPoint);

  Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> extrapolator(
      Acts::EigenStepper<>(m_cfg.magneticField),
      Acts::Navigator({m_cfg.trackingGeometry},
                      logger().cloneWithSuffix("Navigator")),
      logger().cloneWithSuffix("Propagator"));

  // Set up EigenStepper
  Acts::EigenStepper<> stepper(m_cfg.magneticField);

  Acts::PropagatorOptions<Acts::ActionList<Acts::MaterialInteractor>,
                          Acts::AbortList<Acts::EndOfWorldReached>>
      extrapolationOptions(ctx.geoContext, ctx.magFieldContext);
  Acts::PropagatorOptions<> pOptions(ctx.geoContext, ctx.magFieldContext);

  Acts::FullBilloirVertexFitter::Config vertexFitterCfg;
  vertexFitterCfg.extractParameters
      .connect<&Acts::InputTrack::extractParameters>();

  TrackContainer mergedTracks{std::make_shared<Acts::VectorTrackContainer>(),
                              std::make_shared<Acts::VectorMultiTrajectory>()};

  bool first = true;

  auto it1 = m_inputTrackParameters.begin();
  auto it2 = m_inputTrackContainer.begin();

  const auto& itrkcon_tst = *it2;
  const auto& trackContainter_tst = (*itrkcon_tst)(ctx);
  std::uint32_t tipindexOffset = 0;
  for (; it1 != m_inputTrackParameters.end() &&
         it2 != m_inputTrackContainer.end();
       ++it1, ++it2) {
    const auto& itrkpar = *it1;
    const auto& itrkcon = *it2;

    const auto& paramsinput = (*itrkpar)(ctx);
    const auto& trackContainter = (*itrkcon)(ctx);

    if (first) {
      mergedTracks.ensureDynamicColumns(trackContainter);
      first = false;
    }

    auto inputTracks = makeInputTracks(paramsinput);
    int index = -1;

    for (auto& track : inputTracks) {
      index += 1;

      Acts::BoundTrackParameters params1 =
          vertexFitterCfg.extractParameters(track);

      const auto res =
          extrapolator.propagateToSurface(params1, *pSurface, pOptions);

      if (!res.ok()) {
        continue;
      }

      const auto& endParams = *res;

      auto destProxy = mergedTracks.makeTrack();
      auto srcProxy = trackContainter.getTrack(index);

      destProxy.copyFrom(srcProxy, true);
      //destProxy.tipIndex() = srcProxy.tipIndex();//+tipindexOffset;
      destProxy.parameters() = endParams.parameters();
      if (endParams.covariance().has_value()) {
        destProxy.covariance() = endParams.covariance().value();
      }
    }
    tipindexOffset = trackContainter.getTrack(index).tipIndex();
  }

  ActsExamples::ConstTrackContainer outputTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(mergedTracks.container())),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(mergedTracks.trackStateContainer()))};
  
  /*
  ActsExamples::ConstTrackContainer outputTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(mergedTracks.container())),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(trackContainter_tst.trackStateContainer()))};

  */
  m_outputTracks(ctx, std::move(outputTracks));
  return ActsExamples::ProcessCode::SUCCESS;
}