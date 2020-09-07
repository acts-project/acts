// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <memory>

namespace ActsExamples {

class AdaptiveMultiVertexFinderAlgorithm : public ActsExamples::BareAlgorithm {
 public:
  struct Config {
    /// Input track collection
    std::string trackCollection;
  };

  /// Constructor
  AdaptiveMultiVertexFinderAlgorithm(
      const Config& cfg, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Framework execute method
  /// @param [in] context is the Algorithm context for event consistency
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& context) const final override;

 private:
  /// The config class
  Config m_cfg;

  std::vector<Acts::BoundTrackParameters> getInputTrackCollection(
      const ActsExamples::AlgorithmContext& ctx) const;
};

}  // namespace ActsExamples
