// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Vertexing/AdaptiveGridDensityVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts {
class MagneticFieldProvider;
}  // namespace Acts

namespace ActsExamples {
struct AlgorithmContext;

class AdaptiveMultiVertexFinderAlgorithm final : public IAlgorithm {
 public:
  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
  using IPEstimator =
      Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  using Fitter =
      Acts::AdaptiveMultiVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  using Options = Acts::VertexingOptions<Acts::BoundTrackParameters>;

  using VertexCollection =
      std::vector<Acts::Vertex<Acts::BoundTrackParameters>>;

  enum class SeedFinder { GaussianSeeder, AdaptiveGridSeeder };

  struct Config {
    /// Optional. Input track parameters collection
    std::string inputTrackParameters;
    /// Optional. Input trajectories container.
    std::string inputTrajectories;
    /// Output proto vertex collection
    std::string outputProtoVertices;
    /// Output vertex collection
    std::string outputVertices = "vertices";
    /// Enum member determining the choice of the vertex seed finder
    SeedFinder seedFinder;
    /// The magnetic field
    std::shared_ptr<Acts::MagneticFieldProvider> bField;
  };

  AdaptiveMultiVertexFinderAlgorithm(const Config& config,
                                     Acts::Logging::Level level);

  /// Set up vertex seeder and call the function executeAfterSeederChoice.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Find vertices using the adaptive multi vertex finder algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @param seedFinder is the vertex seed finder
  /// @return a process code indication success or failure
  template <typename vseeder_t, typename vfinder_t>
  ProcessCode executeAfterSeederChoice(
      const ActsExamples::AlgorithmContext& ctx,
      const vseeder_t& seedFinder) const;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<std::vector<Acts::BoundTrackParameters>>
      m_inputTrackParameters{this, "InputTrackParameters"};

  ReadDataHandle<TrajectoriesContainer> m_inputTrajectories{
      this, "InputTrajectories"};

  WriteDataHandle<ProtoVertexContainer> m_outputProtoVertices{
      this, "OutputProtoVertices"};

  WriteDataHandle<VertexCollection> m_outputVertices{this, "OutputVertices"};
};
}  // namespace ActsExamples
