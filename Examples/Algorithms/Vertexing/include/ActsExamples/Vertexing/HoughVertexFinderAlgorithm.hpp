// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/EventData/Vertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>

namespace ActsExamples {

class HoughVertexFinderAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Optional. Input space points container.
    std::string inputSpacePoints;
    /// Output vertex collection
    std::string outputVertices;

    /// Configuration for the HoughVertexFinder
    /// See HoughVertexFinder.hpp for description
    std::uint32_t targetSPs = 10000;
    double minAbsEta = 0.3;
    double maxAbsEta = 4.0;
    std::uint32_t minHits = 4;
    Acts::Vector3 defVtxPosition{0., 0., 0.};
  };

  HoughVertexFinderAlgorithm(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// @brief Find a vertex using space points
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};
  WriteDataHandle<VertexContainer> m_outputVertices{this,
                                                    "OutputHoughVertices"};
};

}  // namespace ActsExamples
