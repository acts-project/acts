// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/HoughVertexFinderAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/HoughVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <chrono>
#include <vector>

namespace ActsExamples {

HoughVertexFinderAlgorithm::HoughVertexFinderAlgorithm(const Config& cfg,
                                                       Acts::Logging::Level lvl)
    : IAlgorithm("HoughVertexFinder", lvl), m_cfg(cfg) {
  if (m_cfg.inputSpacepoints.empty()) {
    ACTS_ERROR("You have to provide seeds");
  }
  if (m_cfg.outputVertices.empty()) {
    ACTS_ERROR("Missing output vertices collection");
  }

  m_inputSpacepoints.initialize(m_cfg.inputSpacepoints);
  m_outputVertices.initialize(m_cfg.outputVertices);
}

ProcessCode HoughVertexFinderAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // retrieve input seeds
  const std::vector<SimSpacePoint>& inputSpacepoints = m_inputSpacepoints(ctx);

  Acts::HoughVertexFinder<SimSpacePoint>::Config houghVtxCfg;
  houghVtxCfg.targetSPs = m_cfg.targetSPs;
  houghVtxCfg.minAbsEta = m_cfg.minAbsEta;
  houghVtxCfg.maxAbsEta = m_cfg.maxAbsEta;
  houghVtxCfg.minHits = m_cfg.minHits;
  houghVtxCfg.defVtxPosition = m_cfg.defVtxPosition;
  Acts::HoughVertexFinder<SimSpacePoint> houghVertexFinder(houghVtxCfg);

  // find vertices and measure elapsed time
  auto t1 = std::chrono::high_resolution_clock::now();
  auto vtx = houghVertexFinder.find(inputSpacepoints);
  auto t2 = std::chrono::high_resolution_clock::now();
  if (vtx.ok()) {
    ACTS_INFO("Found a vertex in the event in " << (t2 - t1).count() / 1e6
                                                << " ms");
    ACTS_INFO("Found vertex at x = " << vtx.value()[0]
                                     << "mm, y = " << vtx.value()[1]
                                     << "mm, z = " << vtx.value()[2] << "mm");

    std::vector<Acts::Vertex> vertexCollection;
    vertexCollection.emplace_back(vtx.value());

    // store found vertices
    m_outputVertices(ctx, std::move(vertexCollection));
  } else {
    ACTS_INFO("Not found a vertex in the event after "
              << (t2 - t1).count() / 1e6 << " ms");
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
