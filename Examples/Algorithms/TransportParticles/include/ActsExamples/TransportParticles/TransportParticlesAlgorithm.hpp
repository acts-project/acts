// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>


namespace ActsExamples {
struct AlgorithmContext;

class TransportParticles final : public IAlgorithm {
 public:
  struct Config {
    /// Optional. Input track parameters collection
    std::vector<std::string> inputTrackParameters;
    /// Optional. Input tracks collection
    std::vector<std::string> inputTrackContainer;
    /// Optional. Output tracks collection
    std::string outputTracks;
    /// Output vertex collection
    std::string inputVertices = "vertices";

    bool useRecVtx = false;
    double px = 0;
    double py = 0;
    double pz = 0;

    /// The tracking geometry that should be used.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// The magnetic field that should be used.
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
    /// Extrapolation strategy
    Acts::TrackExtrapolationStrategy extrapolationStrategy =
        Acts::TrackExtrapolationStrategy::firstOrLast;
  };

  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  TransportParticles(Config cfg, Acts::Logging::Level lvl);

  /// Run the ambiguity resolution algorithm.
  ///
  /// @param cxt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

  inline std::vector<Acts::InputTrack> makeInputTracks(
      const TrackParametersContainer& trackParameters) const {
    std::vector<Acts::InputTrack> inputTracks;
    inputTracks.reserve(trackParameters.size());

    for (const auto& trackParam : trackParameters) {
      inputTracks.emplace_back(&trackParam);
    }
    return inputTracks;
  }

 private:
  Config m_cfg;

  std::vector<std::unique_ptr<ReadDataHandle<TrackParametersContainer>>>
      m_inputTrackParameters{};
      
  std::vector<std::unique_ptr<ReadDataHandle<ConstTrackContainer>>>
      m_inputTrackContainer{};

  ReadDataHandle<std::vector<Acts::Vertex>> m_inputVertices{this, "InputVertices"};

  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};


};

}  // namespace ActsExamples
