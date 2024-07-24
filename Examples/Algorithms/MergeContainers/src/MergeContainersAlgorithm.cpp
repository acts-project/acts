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
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <unordered_map>

ActsExamples::MergeContainersAlgorithm::
    MergeContainersAlgorithm(
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
            this,
            "inputTrackParameters#" + std::to_string(m_inputTrackParameters.size())));
    handle->initialize(spName);
  }

  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
}

ActsExamples::ProcessCode
ActsExamples::MergeContainersAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  TrackParametersContainer outputTrackParameter;

  int trk_size = 0;
  for (const auto& itrkpar : m_inputTrackParameters) {
    const auto& handle = (*itrkpar)(ctx);
    trk_size += handle.size();
    std::cout<<"trk_size: "<<handle.size()<<std::endl;
  }
  outputTrackParameter.reserve(trk_size);

  for (const auto& itrkpar : m_inputTrackParameters) {
    const auto& handle = (*itrkpar)(ctx);
    std::copy(handle.begin(), handle.end(), std::back_inserter(outputTrackParameter));
  }
  std::cout<<"merged trk_size: "<<outputTrackParameter.size()<<std::endl;

  m_outputTrackParameters(ctx, std::move(outputTrackParameter));
  return ActsExamples::ProcessCode::SUCCESS;

}