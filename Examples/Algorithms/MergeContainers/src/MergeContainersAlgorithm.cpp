// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/MergeContainers/MergeContainersAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
//#include "Acts/EventData/TrackProxy.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

ActsExamples::MergeContainersAlgorithm::MergeContainersAlgorithm(
    ActsExamples::MergeContainersAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("MergeContainersAlgorithm", lvl),
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

  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTracks) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTracks.emplace_back(
        std::make_unique<ReadDataHandle<ConstTrackContainer>>(
            this, "inputTracks#" + std::to_string(m_inputTracks.size())));
    handle->initialize(spName);
  }

  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode ActsExamples::MergeContainersAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  TrackParametersContainer outputTrackParameters;
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  // std::shared_ptr<Acts::ConstVectorMultiTrajectory> trackStateContainer =
  // m_inputTracks[0](ctx).trackStateContainerHolder();

  TrackContainer mergedTracks{trackContainer, trackStateContainer};

  int trk_size = 0;
  for (const auto& itrkpar : m_inputTrackParameters) {
    const auto& handle = (*itrkpar)(ctx);
    trk_size += handle.size();
    // std::cout << "trk_size: " << handle.size() << std::endl;
  }
  outputTrackParameters.reserve(trk_size);

  for (const auto& itrkpar : m_inputTrackParameters) {
    const auto& handle = (*itrkpar)(ctx);
    std::copy(handle.begin(), handle.end(),
              std::back_inserter(outputTrackParameters));
  }

  using TrackProxyType = Acts::TrackContainer<Acts::ConstVectorTrackContainer,
                                              Acts::ConstVectorMultiTrajectory,
                                              std::shared_ptr>::ConstTrackProxy;

  trk_size = 0;

  for (const auto& itrkpar : m_inputTracks) {
    const auto& handle = (*itrkpar)(ctx);
    mergedTracks.ensureDynamicColumns(handle);
    break;
    //trk_size += handle.size();
  }

  //mergedTracks.reserve(trk_size);

  uint32_t tipIndex = 0;
  for (const auto& itrkpar : m_inputTracks) {

    uint32_t tipIndexTmp = 0;
    const auto& handle = (*itrkpar)(ctx);
    for (auto itrack = handle.begin(); itrack != handle.end(); ++itrack) {
      auto destProxy = mergedTracks.makeTrack();
      destProxy.copyFrom((*itrack), true);
      destProxy.tipIndex() = (*itrack).tipIndex();
      destProxy.stemIndex() = (*itrack).stemIndex();
      tipIndexTmp = (*itrack).tipIndex();
    }
    tipIndex += tipIndexTmp;
  }

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};

  m_outputTrackParameters(ctx, std::move(outputTrackParameters));
  m_outputTracks(ctx, std::move(constTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}